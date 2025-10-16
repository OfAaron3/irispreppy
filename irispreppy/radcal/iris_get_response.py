import pickle
import urllib.error
import urllib.request
from datetime import datetime as dt
from glob import glob
from os import path, remove

import numpy as np
from astropy.time import Time
from bs4 import BeautifulSoup
from scipy.interpolate import interp1d
from scipy.io import readsav



def iris_get_response(date=None, version=0, response_file=None, pre_launch=False, full=False, angstrom=False, update_response=True, quiet=False):
    """
    Functions identically to that of `iris_get_response.pro`. Generates the IRIS response required for radiometric calibration.

    Parameters
    ----------
    date : string
        Date of observation. e.g. '2013-06-28T02:27:46.00Z' Default is now.
    version : int
        Which version of the response file you want. Default is newest. Version numbers are 1 indexed.  Default: 0 (i.e. the newest).
    response_file : string
        Name of the response file you want to use. Must exactly match, including extension. Else reverts to default response of the most recent.
    pre_launch : bool 
        Analaguous to version=2. Default: False.
    full : bool
        Full effective area structure is returned with cryptic coefficients. Default: False.
    angstrom : bool
        Whether to return lambda in angstroms. If False, lambda is returned in nm. Default: False.
    update_response: bool
        Whether to call `update_response_files(quiet=True)` to check for new response files before starting calibration. Default: True.
    quiet : bool
        Whether to suppress all print statements including messages about failing to find specified response file and reverting to default. Default: False.

    Returns
    -------

    o : dict
        IRIS response. A dict was chosen to mimic an IDL struct.

    Notes
    -----
        In order to mimic `iris_get_response.pro`, there are three different parameters that specify the version of the response file to use. If none are set, the most recent response file will be used. Their precedence (if all are set) is:
            `pre_launch` > `version` > `response_file`
    
    Example
    -------
    >>> from astropy.io import fits
    >>> import irispreppy as ip
    >>> f=fits.open('iris_raster.fits')
    >>> frc=ip.iris_get_response(f[0].header['DATE_OBS'])
    """

    if date is None:
        date=dt.strftime(dt.now(), '%Y-%m-%dT%H:%M:%S.%fZ')
    toppath=path.dirname(path.realpath(__file__))
    resppath=path.join(toppath, "responses")
    resps=glob(path.join(resppath, "*.pkl"))
    
    if update_response:
        new=update_response_files(quiet=True)
        if new:
            resps=glob(path.join(resppath, "*.pkl"))
    
    resps.sort()

    #0a Opening correct file
    if pre_launch:
        response=resps[1] #This is version 2, which is file 1.
    elif version!=0:
        if version<=0:
            if not quiet:
                print("No such version of response file. Defaulting to most recent version.")
            response=resps[-1]
        elif version<=len(resps):
            response=resps[version-1]
        else:
            if not quiet:
                print("Requested version of response file not found. Defaulting to most recent version.")
            response=resps[-1]
    elif response_file is not None:
        if toppath+"/responses/"+response_file not in resps:
            if not quiet:
                print(response_file+" not found. Using most recent file.")
            response=resps[-1]
        else:
            response="./responses/"+response_file
    else:
        response=resps[-1]

    with open(response, "rb") as pklin:
        r=pickle.load(pklin) #Loading in the response file and calling it r

    #0b Handling keywords
    if int(r['version'])<2:
        if angstrom:
            r['lambda']=r['lambda']*10
        return(r)

    ntt=1 #ntt is always 1 since this will have a wrapper around it
    
    #1. Output structure
    if full:
        o1=r #Output 1. Is temporary
    else:
        keys=['LAMBDA', 'AREA_SG', 'NAME_SG', 'DN2PHOT_SG', 'AREA_SJI', 'NAME_SJI', 'DN2PHOT_SJI', 'COMMENT', 'VERSION', 'VERSION_DATE']
        o1={k:r[k] for k in keys}
        o1['DATE-OBS']=date

    o1['AREA_SG']=np.zeros_like(o1['AREA_SG'])
    o1['AREA_SJI']=np.zeros_like(o1['AREA_SJI'])
    o=[o1 for i in range(0, ntt)] #Output
    del o1

    #2. FUV SG Effective Areas
    # ; Rough SG spectral ranges. Setting eff.area to 0 outside of these
    lamran=[[133.1,135.9],[138.8,140.8]]
    #Not entirley sure how the coeff array is organised, but the comment on the IDL version says,
    # "; time dependent response for sz[3]=3 wavelengths". The index in Python is [0] though
    sz=r['coeffs_fuv'].shape 
    rr=np.zeros((ntt, sz[0]))
    for j in range(0, sz[0]):
        rr[:,j]=fit_iris_xput_lite(date, r['c_f_time'], r['coeffs_fuv'][j])
    #; interpolate onto lambda grid, separately for each of the two FUV CCDs
    for j in range(0, 2):
        w=(r['lambda']>=lamran[j][0]) & (r['lambda']<=lamran[j][1])
        for k in range(0, ntt):
            interp=interp1d(r['c_f_lambda'][j:j+2], rr[k, j:j+2], fill_value='extrapolate')
            #If you feel uneasy about this extrapolation, this is how iris_get_resposne.pro works implicitly
            o[k]['AREA_SG'][0,w]=interp(r['lambda'][w])
            # ;Version 009+ only: Remove wavelength dependence for the Si IV part of the FUV window.
            if int(o['VERSION'])>=9 and j==1:
                o[k]['AREA_SG'][0, w]=np.ones(np.sum(w))*np.mean(o[k]['AREA_SG'][0,w])


    #3. NUV SG Effective Areas
    # ; Rough SG spectral ranges. Setting eff.area to 0 outside of these
    lamran=[278.2,283.5]
    #Not entirley sure how the coeff array is organised, but the comment on the IDL version says,
    # "; time dependent response for sz[3]=3 wavelengths". The index in Python is [0] though
    sz=r['coeffs_nuv'].shape 
    rr=np.zeros((ntt, sz[0]))
    for j in range(0, sz[0]):
        rr[:,j]=fit_iris_xput_lite(date, r['c_n_time'], r['coeffs_nuv'][j])

    # ; apply wavelength-independent factor to all wavelengths sz[3]
    if int(o1['VERSION'])==7:
        #;determine if input time contains period of A1 QS 2820-2832A trend...
        #This appears to be a "quick fix" for the effective area drop
        trend_tim=np.array([t[:10] for t in r['TREND_TIM'].astype(str)])
        if (trend_tim==date[:10]).any():
            tt1=np.argmax(trend_tim==date[:10])
            if (tt1==0 and trend_tim[0]==date[:10]) or tt1>0:
                trendy=r['TREND'][tt1]
                for j in range(0, sz[0]):
                    rr[:,j]=rr[:,j]*trendy                     

    #; interpolate onto lambda grid
    w=(r['lambda']>=lamran[0]) & (r['lambda']<=lamran[1])
    if int(r['version'])<3:
        for k in range(0, ntt):
            interp=interp1d(r['c_n_lambda'], rr[k], fill_value='extrapolate')
            o[k]['AREA_SG'][1,w]=interp(r['lambda'][w])
    else: #I guess for version>=3, len(r['c_n_lambda'])>2
        interp=interp1d(r['c_n_lambda'], rr[k], fill_value='extrapolate', kind='quadratic')
        o[k]['AREA_SG'][1,w]=interp(r['lambda'][w])


    #4. SJI Effective Areas
    if int(r['version'])<3:
        sz=r['coeffs_sji'].shape 
        for j in range(0, sz[0]):
        # ; calculate pre-launch area from the individual elements
            pl_a=r['geom_area']
            for k in range(0, len(r['index_el_sji'][0])):
                pl_a=(pl_a*np.array([[r['elements'][r['index_el_sji'][0]][i][1]] for i in range(0, r['elements'][r['index_el_sji'][0]].shape[0])])).T
                pl_a=pl_a[:,0,:] #Because of the way I get this to work, I introduce an extra 1-length axis
            rr=fit_iris_xput_lite(date, r['c_s_time'][j], r['coeffs_sji'][j])
            for k in range(0, ntt):
                o[k]['AREA_SJI'][j]=pl_a*rr[k]

    else:
        for nuv in range(0, 2):
        # ; calculate baseline SJI area curves
            asji=r['geom_area']
            for k in range(0, len(r['index_el_sji'][nuv*2])):
                arr=np.array([[r['elements'][r['index_el_sji'][2:4, 3]][i]] for i in range(0, r['elements'][r['index_el_sji'][2:4, 3]].shape[0])])
                asji=asji*(np.array([arr[0][0][1], arr[1][0][1]])).T
                del arr
            # ; apply time dependent profile shape adjustment to FUV SJI
            if ~nuv:
                # ; FUV: apply FUV SG "slant", then normalize so that a weighted (2.4:1)
                # ;      sum at C II and Si IV gives constant response
                wei=[2.4,1] # ; typical solar ratio CII : SiIV
                wav=r['c_f_lambda']
                nwv=len(wav)
                wav=[wav[0], (wav[nwv-2]*2+wav[nwv-1])/3] # ; 2 wvlngts in nm
                # ; calculate baseline SG area for scaling purposes
                asg=r['geom_area']
                for k in range(0, len(r['index_el_sg'][nuv])):
                    asg=asg*r['elements'][r['index_el_sg']][nuv, k][1].T
                # ; SG and SJI areas at wav
                interp=interp1d(r['lambda'], asg, fill_value='extrapolate')
                asg2=interp(wav)
                asj2=np.zeros((2,2))
                for j in range(0,2):
                    interp=interp1d(r['lambda'], asji[:,j], fill_value='extrapolate')
                    asj2[:,j]=interp(wav)
                # ; calculate the normalized slant function scal, apply to asji
                for k in range(0, ntt):
                    # ; ; best-estimate slant, i.e., eff.area @ wav / baseline SG @ wav
                    interp=interp1d(o[k]['LAMBDA'], o[k]['AREA_SG'][0], fill_value='extrapolate')
                    sca2=interp(wav)/asg2
                    # ; normalize slant so that total(wei*asj2*sca2)/total(wei*asj2)=1
                    for j in range(0, 2):
                        sca2n=sca2*np.sum(wei*asj2[:,j], axis=None)/np.sum(wei*asj2[:,j]*sca2)
                        interp=interp1d(wav, sca2n, fill_value='extrapolate')
                        scaln=interp(r['lambda'])
                        o[k]['AREA_SJI'][j]=asji[:,j]*scaln
            else:
                # ; NUV: essentially same calculation as r.version=3
                for k in range(0, ntt):
                    o[k]['AREA_SJI'][2:4]=asji
            
        for j in range(0, 4):
            # ; SJI specific time dependency
            rr=fit_iris_xput_lite(date, r['c_s_time'][j], r['coeffs_sji'][j])
            for k in range(0, ntt):
                o[k]['AREA_SJI'][j]=o[k]['AREA_SJI'][j]*rr[k]

    if angstrom:
        o['lambda']=o['lambda']*10
    return(o)


def fit_iris_xput_lite(tt0, tcc0, ccc):
    '''
    Stripped down form of `fit_iris_xput.pro`, with only the sections relevant to `iris_get_response.pro`.
    I am so sorry, but I have no idea what the last two inputs are.

    Parameters
    ----------
    tt0 : string
        Observation date.
    tcc0 : array_like
        Time Coefficients.
    ccc :  array_like
        Coefficients.
        
    Returns
    -------
    f : float
       IRIS crossput.

    '''
    tex=1.5 # ; exponent for transition between exp.decay intervals
    if tcc0.shape[1]!=2:
        raise RuntimeError("Incorrect number of elements in tcoef (tcc0)")
    m=tcc0.shape[0]
    #This is crazy. Originally here they did 
    #m=size(tcc0); if m[1] ne 2; return, 0; endif; m=m[m[0]+2]/2.

    #; times in years since tcc[0,0]
    #The original is using NASA's epoch of 1979. I'm using 1970, as is standard.
    nasaLag=Time('1979-01-01T00:00:00.00Z').unix
    yr2sec=365.25*24*60*60
    tcc=(tcc0+nasaLag)/yr2sec 
    tt=Time(tt0).unix/yr2sec
    ntt=1 
    #Always going to be 1 in this instance. Script can, in theory, take more than one input.
         
    if ccc.size!=3*m:
        raise RuntimeError("Incorrect number of elements in tcoef (tcco)")
    # ; calculation of output fit
    ee=ccc[:,2] 
    a=np.zeros((ntt, 2*m)) #; base vector
    #no need to reform in python, it's 2D by default here
    for j in range(0, m):
        # ; base vector for interval of constant multiplier
        if j>0:
            nww=int((tt>=tcc[j, 0]) & (tt<=tcc[j, 1]))
            ww=nww-1
        else:
            nww=int(tt<=tcc[j,1])
            ww=nww-1

        if nww>0:
            a[ww, 2*j]=1
            a[ww, 2*j+1]=np.exp(ee[j]*(tt-tcc[j,0]))
        # ; base vector for interval when multiplier is linear function of time
        # Sometimes dtt<0, so have to NaN it before the power to stop a warning
        if j>0:
            nww=int((tt>tcc[j-1,1]) & (tt < tcc[j,0]))
            ww=nww-1
            if nww>0:
                dtt=(tt-tcc[j-1,1])/(tcc[j,0]-tcc[j-1, 1])
                if dtt<0:
                    dtt=np.nan
                a[ww, 2*j]=dtt**tex
                a[ww, 2*j+1]=dtt**tex*np.exp(ee[j]*(tt-tcc[j,0]))
        if j < (m-1):
            nww=int((tt>tcc[j,1]) & (tt<tcc[j+1, 0]))
            ww=nww-1
            if nww>0:
                dtt=(tt-tcc[j,1])/(tcc[j+1,0]-tcc[j,1])
                if dtt<0:
                    dtt=np.nan
                a[ww, 2*j]=1-dtt**tex
                a[ww, 2*j+1]=(1-dtt**tex)*np.exp(ee[j]*(tt-tcc[j,0]))

    cc=ccc[:,:2].flatten()
    f=a@cc
    return(f)


def update_response_files(quiet=False):
    '''
    Checks https://hesperia.gsfc.nasa.gov/ssw/iris/response/ for new IRIS response files, and downloads any it finds.

    Parameters
    ----------
    quiet : bool
        Whether to suppress printing status, and URL errors. Default: False.
    
    Returns
    -------
    new : bool
        Whether any new response files were found.

    Example
    -------
    >>> import irispreppy as ip
    >>> f=ip.update_response_files()
    '''
    toppath=path.dirname(path.realpath(__file__))
    resppath=path.join(toppath, "responses")
    resps=glob(path.join(resppath, "*.pkl"))
    resps.sort()
    new=False
    try:
        if not quiet:
            print("Connecting to hesperia...")
        with urllib.request.urlopen("https://hesperia.gsfc.nasa.gov/ssw/iris/response/") as respurl:
            htmlsoup=BeautifulSoup(respurl, 'html.parser')
        for tags in htmlsoup.find_all('a'):
            href=tags.get('href')
            if "sra" in href and path.join(resppath, href[:-4]+'pkl') not in resps:
                new=True
                if not quiet:
                    print("New response file found, "+href+'.\nDownloading...')
                urllib.request.urlretrieve("https://hesperia.gsfc.nasa.gov/ssw/iris/response/"+href, "temp.geny")
                newgeny=readsav('temp.geny')
                remove('temp.geny')
                recgeny=newgeny[list(newgeny.keys())[0]][0]
                with open(toppath+"/responses/"+href[:-4]+'pkl', "wb") as pklout:
                    pickle.dump(recgeny, pklout)
                
    except urllib.error.URLError:
        if not quiet:
            print("Either hesperia is unreachable or you are not connected to the internet. Cannot check for new response files.")
    except:
        if not quiet:
            print("Hesperia is reachable but not behaving expectedly. Cannot check for new response files.")

    if not quiet and not new:
            print("No new response files found.")

    return(new)
