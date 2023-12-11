import pickle
import urllib.error
import urllib.request
from datetime import datetime as dt
from glob import glob as ls
from os import path, remove

import numpy as np
from astropy.time import Time
from bs4 import BeautifulSoup
from scipy.interpolate import interp1d
from scipy.io import readsav


def iris_get_response(date=dt.strftime(dt.now(), '%Y-%m-%dT%H:%M:%S.%fZ'), version=0, response_file=None, pre_launch=False, full=False, angstrom=False, quiet=False):
    '''Intended to use in place of iris_get_response.pro
    Input Parameters:
        date: Time or time list. Default is now.
        version: Which version of the response file you want. Default is newest. Version numbers are 1 indexed, so the default 0 becomes -1.
        response_file: Name of the response file you want to use. Must exactly match, including extension.
        pre_launch: Not sure why this is in the original, but it is analaguous to version=2. Default is False.
        full: Full effective area structure is returned with cryptic coefficients. Default is False.
        angstrom: If True, lambda is returned in angstroms. If False, lambda is retunred in nm. Default is False.
        quiet: If true, prints messages about contacting hesperia

    Notes:
        1. version, response_file, and prelaunch all perform the same function, here is their precedence,
        pre_launch>version>response_file 
        2. Code automatically checks https://hesperia.gsfc.nasa.gov/ssw/iris/response/ for new response files. If this url
        changes in the future, do a search and replace. The files are assumed to be geny IDL structs.
        3. All original comments will be preceeded by ;, as is convention in IDL
        4. Translated from iris_get_response.pro. Originally by J.P.Weulser, S.L. Freeland, and G.Chintzoglou

    History:
        2021-12-14 - A.W.Peat - Translated from IDL and added QOL improvements
    '''

    toppath=path.dirname(path.realpath(__file__))
    resppath=path.join(toppath, "responses")
    resps=ls(path.join(resppath, "*.pkl"))
    resps.sort()

    #Checks for new responses everytime it is run
    try:
        with urllib.request.urlopen("https://hesperia.gsfc.nasa.gov/ssw/iris/response/") as respurl:
            htmlsoup=BeautifulSoup(respurl, 'html.parser')
        for tags in htmlsoup.find_all('a'):
            href=tags.get('href')
            if "sra" in href and path.join(resppath, href[:-4]+'pkl') not in resps:
                if not quiet:
                    print("New response file found, "+href+'.\nDownloading...')
                urllib.request.urlretrieve("https://hesperia.gsfc.nasa.gov/ssw/iris/response/"+href, "temp.geny")
                newgeny=readsav('temp.geny')
                remove('temp.geny')
                recgeny=newgeny[list(newgeny.keys())[0]][0]
                with open(toppath+"/responses/"+href[:-4]+'pkl', "wb") as pklout:
                    pickle.dump(recgeny, pklout)

                resps=ls(toppath+"/responses/*.*") #Needs to reload responses if a new one is found
                resps.sort()
    except urllib.error.URLError:
        if not quiet:
            print("You are not connected to the internet. Cannot check for new response files.")
    except:
        if not quiet:
            print("Hesperia is reachable but not loading. Cannot check for new response files.")

    #0a Opening correct file
    if pre_launch:
        response=resps[1] #0 indexing
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
    elif response_file!=None:
        if toppath+"/responses/"+response_file not in resps:
            if not quiet:
                print(response_file+" not found. Using most recent file.")
            response=resps[-1]
        else:
            response="./responses/"+response_file
    else:
        response=resps[version-1]

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
        w=np.where((r['lambda']>=lamran[j][0]) & (r['lambda']<=lamran[j][1]))
        for k in range(0, ntt):
            interp=interp1d(r['c_f_lambda'][j:j+2], rr[k, j:j+2], fill_value='extrapolate')
            #If you feel uneasy about this extrapolation, this is how iris_get_resposne.pro works implicitly
            o[k]['AREA_SG'][0,w]=interp(r['lambda'][w])


    #3. NUV SG Effective Areas
    # ; Rough SG spectral ranges. Setting eff.area to 0 outside of these
    lamran=[278.2,283.5]
    #Not entirley sure how the coeff array is organised, but the comment on the IDL version says,
    # "; time dependent response for sz[3]=3 wavelengths". The index in Python is [0] though
    sz=r['coeffs_nuv'].shape 
    rr=np.zeros((ntt, sz[0]))
    for j in range(0, sz[0]):
        rr[:,j]=fit_iris_xput_lite(date, r['c_n_time'], r['coeffs_nuv'][j])
    #; interpolate onto lambda grid
    w=np.where((r['lambda']>=lamran[0]) & (r['lambda']<=lamran[1]))
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
            rr=fit_iris_xput_lite(time, r['c_s_time'][j], r['coeffs_sji'][j])
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
    Stripped down form of fit_iris_xput.pro, using only the things 
    get_iris_response.pro uses.
    I am so sorry, but I have no idea what any of these keywords are.
    The previous documentation is very cryptic. I will include ALL of their comments.

    Notes:
        1. All original comments will be preceeded by ;, as is convention in IDL
        2. Based on fit_iris_xput.pro by JPW.

    History:
        2021-12-14 - A.W.Peat - Translated from IDL
    '''
    tex=1.5 # ; exponent for transition between exp.decay intervals
    if tcc0.shape[1]!=2:
        raise RuntimeError("Incorrect number of elements in tcoef (tcco)")
    m=tcc0.size//2
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
            ww=np.where((tt>=tcc[j, 0]) & (tt<=tcc[j, 1]))
            nww=len(ww)
        else:
            ww=np.where(tt<=tcc[j,1])
            nww=len(ww)

        if nww>0:
            a[ww, 2*j]=1
            a[ww, 2*j+1]=np.exp(ee[j]*(tt-tcc[j,0]))
        # ; base vector for interval when multiplier is linear function of time
        # Sometimes dtt<0, so have to NaN it before the power to stop a warning
        if j>0:
            ww=np.where((tt>tcc[j-1,1]) & (tt < tcc[j,0]))
            nww=len(ww)
            if nww>0:
                dtt=(tt-tcc[j-1,1])/(tcc[j,0]-tcc[j-1, 1])
                if dtt<0:
                    dtt=np.nan
                a[ww, 2*j]=dtt**tex
                a[ww, 2*j+1]=dtt**tex*np.exp(ee[j]*(tt-tcc[j,0]))
        if j < (m-1):
            ww=np.where((tt>tcc[j,1]) and (tt<tcc[j+1, 0]))
            nww=len(ww)
            if nww>0:
                dtt=(tt-tcc[j,1])/(tcc[j+1,0]-tcc[j,1])
                if dtt<0:
                    dtt=np.nan
                a[ww, 2*j]=1-dtt**tex
                a[ww, 2*j+1]=(1-dtt**tex)*np.exp(ee[j]*(tt-tcc[j,0]))

    cc=ccc[:,:2].flatten()
    f=a@cc
    return(f)


if __name__=="__main__":
    #When script is called directly, it just looks for new response files#
    toppath=path.dirname(path.realpath(__file__))
    resppath=path.join(toppath, "responses")
    resps=ls(path.join(resppath, "*.pkl"))
    resps.sort()
    new=False
    try:
        with urllib.request.urlopen("https://hesperia.gsfc.nasa.gov/ssw/iris/response/") as respurl:
            htmlsoup=BeautifulSoup(respurl, 'html.parser')
        for tags in htmlsoup.find_all('a'):
            href=tags.get('href')
            if "sra" in href and path.join(resppath, href[:-4]+'pkl') not in resps:
                print("New response file found, "+href+'.\nDownloading...')
                urllib.request.urlretrieve("https://hesperia.gsfc.nasa.gov/ssw/iris/response/"+href, "temp.geny")
                newgeny=readsav('temp.geny')
                remove('temp.geny')
                recgeny=newgeny[list(newgeny.keys())[0]][0]
                with open(toppath+"/responses/"+href[:-4]+'pkl', "wb") as pklout:
                    pickle.dump(recgeny, pklout)

                resps=ls(toppath+"/responses/*.*") #Needs to reload responses if a new one is found
                resps.sort()
    except urllib.error.URLError:
        print("You are not connected to the internet. Cannot check for new response files.")
        new=True
    except:
        print("Hesperia is reachable but not loading. Cannot check for new response files.")
        new=True
    if not new:
        print("No new response files found.")

    
