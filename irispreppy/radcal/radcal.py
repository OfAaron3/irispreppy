import datetime as dt
from copy import deepcopy as dc

from astropy.io import fits
import numpy as np
import scipy.stats as scist
from weno4 import weno4

from . import iris_get_response as igr

#################################################################

def radiometric_calibrate(ras, quiet=True, error=False):
    '''
    Performs radiometric calibration of IRIS spectrograph files.\n
    (See Section 5.2 of ITN26 for more details - https://iris.lmsal.com/itn26/calibration.html#radiometric-calibration)
    
    Parameters
    ----------
    ras : astropy.io.fits.hdu.hdulist.HDUList
        Input IRIS raster
    quiet : bool
        Whether to suppress print statements. Default: True
    error : bool
        Whether to calculate errors (beta). Ignored if `ras` is a full disc mosaic. Default: False

    Returns
    -------
    hdul : astropy.io.fits.hdu.hdulist.HDUList
        The calibrated data
    hdule : astropy.io.fits.hdu.hdulist.HDUList
        If error=True. The poissonian error and dark current on the data (beta)

    Example
    -------
    >>> from astropy.io import fits
    >>> import irispreppy as ip
    >>> f=fits.open('iris_raster.fits')
    >>> frc=ip.radiometric_calibrate(f)

    '''

    if type(ras)==fits.hdu.hdulist.HDUList:
        if 'TELESCOP' not in ras[0].header:
            if not quiet:
                print("Telescope keyword not present. Assuming full disc mosaic.")
        else:
            assert ras[0].header['TELESCOP']=='IRIS'
        rasfits=ras

    else:
        raise ValueError("Must supply astropy.io.fits.hdu.hdulist.HDUList of an IRIS observation")


    ############
    # Response #
    ############

    rasfits.verify('fix')

    if 'STARTOBS' not in rasfits[0].header:
        begin=dt.datetime.strptime(rasfits[0].header['DATE_OBS'], '%Y-%m-%dT%H:%M:%S.%f')
    else:
        begin=dt.datetime.strptime(rasfits[0].header['STARTOBS'], '%Y-%m-%dT%H:%M:%S.%f')
    if 'ENDOBS' not in rasfits[0].header:
        end=dt.datetime.strptime(rasfits[0].header['DATE_END'], '%Y-%m-%dT%H:%M:%S.%f')
    else:
        end=dt.datetime.strptime(rasfits[0].header['ENDOBS'], '%Y-%m-%dT%H:%M:%S.%f')

    midtime=dt.datetime.strftime((begin+((end-begin)/2)), '%Y-%m-%dT%H:%M:%S.%fZ')

    response=(igr.iris_get_response(midtime, quiet=quiet))[0]

    ###################################
    # Find which index is FUV and NUV #
    ###################################

    if response['NAME_SG'][0]==b'FUV' and response['NAME_SG'][1]==b'NUV':
        FUVind=0
        NUVind=1
    elif response['NAME_SG'][1]==b'FUV' and response['NAME_SG'][0]==b'NUV':
        FUVind=1
        NUVind=0
    else:
        print("[NAME_SG]="+str(response['NAME_SG']))
        raise RuntimeError("FUV and NUV cannot be found automatically. Please check ['NAME_SG'] from irisresponse above.")

    ##################################
    # Find where is what in the fits #
    ##################################

    indices={rasfits[0].header[name]: ind+1 for ind, name in enumerate(rasfits[0].header['TDESC*'])}

    if indices=={}:
        #Full disc mosaic
        if rasfits[0].header['CRVAL3']>2000:
            indices={'fdNUV':0}
        else:
            indices={'fdFUV':0}    

    ##########
    # Blanks #
    ##########

    if not list(indices.keys())[0][:2]=='fd':
        #Blanks not to be removed from fd mosiacs
        #Piggybacking off the last if a little, but it's probably okay
        blanks={}
        for key in indices:
            blanks[key]=[]
            for yy in range(0, rasfits[indices[key]].header['NAXIS2']):
                if (rasfits[indices[key]].data[:,yy]==-200).all():
                    blanks[key].append(yy)


    ###################
    # Wavelength axes #
    ###################

    FUV=np.where(response['AREA_SG'][FUVind]>0)[0]
    FUVcutoff=np.where((FUV[1:]-FUV[:-1])>1)[0][0]+1

    FUV1=FUV[:FUVcutoff]
    FUV1=response['LAMBDA'][FUV1]*10

    FUV2=FUV[FUVcutoff:]
    FUV2=response['LAMBDA'][FUV2]*10

    NUV=(response['LAMBDA'][response['AREA_SG'][NUVind]>0])*10

    #################
    # Photon Energy #
    #################
    #Lambda is in A. 1e7 ergs = 1 Joule
    h=6.62607004e-34
    c=3e8
    eFUV1=1e7*h*c/(FUV1*1e-10)
    eFUV2=1e7*h*c/(FUV2*1e-10)
    eNUV=1e7*h*c/(NUV*1e-10)

    ##################
    # Effective Area #
    ##################

    aFUV1=np.trim_zeros(response['AREA_SG'][FUVind][FUV[:FUVcutoff]])
    aFUV2=np.trim_zeros(response['AREA_SG'][FUVind][FUV[FUVcutoff:]])
    aNUV=np.trim_zeros(response['AREA_SG'][NUVind])

    del FUV

    ###########
    # DN2PHOT #
    ###########
    d2pFUV=response['DN2PHOT_SG'][FUVind]
    d2pNUV=response['DN2PHOT_SG'][NUVind]

    ###############
    #Exposure Time#
    ###############
    if 'EXPTIMEN' in rasfits[0].header:
        tnuv=rasfits[0].header['EXPTIMEN']
    else:
        tnuv=rasfits[0].header['EXPTIME']
    if 'EXPTIMEF' in rasfits[0].header:
        tfuv=rasfits[0].header['EXPTIMEF']
    else:
        tfuv=rasfits[0].header['EXPTIME']

    ############
    #Slit Width#
    ############
    wslit=np.pi/(180*3600*3)

    ##################################
    #Spectral Pixel Width [angstroms]#
    ##################################

    if 'fdNUV' in indices or 'fdFUV' in indices:
        pixl={key: rasfits[indices[key]].header['CDELT3'] for key in indices}
    else:
        pixl={key: rasfits[indices[key]].header['CDELT1'] for key in indices}

    ##############################
    #Spatial Pixel Size [radians]#
    ##############################

    #ITN26 is a little cryptic about this, but I swear this is right
    pixxy={key: rasfits[indices[key]].header['CDELT2']*np.pi/(180*3600) for key in indices}

    #############
    # Constants #
    #############

    const={}
    for key in indices:
        if key=='fdFUV':
            const[key]=d2pFUV/(pixxy[key]*pixl[key]*tfuv*wslit)
        elif key=='fdNUV':
                const[key]=d2pNUV/(pixxy[key]*pixl[key]*tnuv*wslit)
        elif 'FUV' in rasfits[0].header['TDET'+str(indices[key])]:
            const[key]=d2pFUV/(pixxy[key]*pixl[key]*tfuv*wslit)
        else:
            const[key]=d2pNUV/(pixxy[key]*pixl[key]*tnuv*wslit)

    ###########################################################
    # Wavelength Trimming and Radiometric Calibration Factors #
    ###########################################################
    
    lamwin={} #LAMbda WINdow
    wvlns={}
    rcfs={} #Radiometric Calibration FactorS
    for key in indices:
        if key=='fdFUV':
            wvlns[key]=(np.arange(0, rasfits[indices[key]].header['NAXIS3'])-rasfits[indices[key]].header['CRPIX3']+1)*rasfits[indices[key]].header['CDELT3']+rasfits[indices[key]].header['CRVAL3']
            lamwin[key]=np.arange(0, len(wvlns[key]))[(wvlns[key]>FUV1[0])&(wvlns[key]<FUV1[-1])]
            lamwin[key]=lamwin[key][0:len(lamwin[key]):len(lamwin[key])-1]            
            rcfs[key]=weno4(wvlns[key][lamwin[key][0]:lamwin[key][1]], FUV1, eFUV1/aFUV1)*const[key]

        elif key=='fdNUV':
            wvlns[key]=(np.arange(0, rasfits[indices[key]].header['NAXIS3'])-rasfits[indices[key]].header['CRPIX3']+1)*rasfits[indices[key]].header['CDELT3']+rasfits[indices[key]].header['CRVAL3']
            lamwin[key]=np.arange(0, len(wvlns[key]))[(wvlns[key]>NUV[0])&(wvlns[key]<NUV[-1])]
            lamwin[key]=lamwin[key][0:len(lamwin[key]):len(lamwin[key])-1]            
            rcfs[key]=weno4(wvlns[key][lamwin[key][0]:lamwin[key][1]], NUV, eNUV/aNUV)*const[key]

        elif rasfits[0].header['TDET'+str(indices[key])]=='FUV1':
            wvlns[key]=(np.arange(0, rasfits[indices[key]].header['NAXIS1'])-rasfits[indices[key]].header['CRPIX1']+1)*rasfits[indices[key]].header['CDELT1']+rasfits[indices[key]].header['CRVAL1']
            lamwin[key]=np.arange(0, len(wvlns[key]))[(wvlns[key]>FUV1[0])&(wvlns[key]<FUV1[-1])]
            lamwin[key]=lamwin[key][0:len(lamwin[key]):len(lamwin[key])-1]            
            rcfs[key]=weno4(wvlns[key][lamwin[key][0]:lamwin[key][1]], FUV1, eFUV1/aFUV1)*const[key]

        elif rasfits[0].header['TDET'+str(indices[key])]=='FUV2':
            wvlns[key]=(np.arange(0, rasfits[indices[key]].header['NAXIS1'])-rasfits[indices[key]].header['CRPIX1']+1)*rasfits[indices[key]].header['CDELT1']+rasfits[indices[key]].header['CRVAL1']
            lamwin[key]=np.arange(0, len(wvlns[key]))[(wvlns[key]>FUV2[0])&(wvlns[key]<FUV2[-1])]
            lamwin[key]=lamwin[key][0:len(lamwin[key]):len(lamwin[key])-1]            
            rcfs[key]=weno4(wvlns[key][lamwin[key][0]:lamwin[key][1]], FUV2, eFUV2/aFUV2)*const[key]

        elif rasfits[0].header['TDET'+str(indices[key])]=='NUV':
            wvlns[key]=(np.arange(0, rasfits[indices[key]].header['NAXIS1'])-rasfits[indices[key]].header['CRPIX1']+1)*rasfits[indices[key]].header['CDELT1']+rasfits[indices[key]].header['CRVAL1']
            lamwin[key]=np.arange(0, len(wvlns[key]))[(wvlns[key]>NUV[0])&(wvlns[key]<NUV[-1])]
            lamwin[key]=lamwin[key][0:len(lamwin[key]):len(lamwin[key])-1]            
            rcfs[key]=weno4(wvlns[key][lamwin[key][0]:lamwin[key][1]], NUV, eNUV/aNUV)*const[key]

        else:
            raise ValueError("You have detectors that are not FUV1, FUV2, or NUV in your fits file.")

    for i in list(rcfs.keys()):
        rcfs[i]=rcfs[i].astype(np.float32)
    #Gotta change them to 32 bit floats. They are 64!!

    ###########################
    # Radiometric Calibration #
    ###########################
    if not quiet:
        print("Creating calibrated spectrograph fits...")    
        print("(This may take some time as it needs to recalculate the excessive amount of stats in the headers)")


    hdr0=dc(rasfits[0].header)
    hdr0['HISTORY']='NUV and FUV radiometric calibration performed on '+dt.datetime.now().strftime("%Y-%m-%d")
    hdr0['HISTORY']='FITS made with astropy on '+dt.datetime.now().strftime("%Y-%m-%d")
    hdr0['BUNIT']="erg s^-1 cm^-2 angstrom^-1 sr^-1"
    dat={}
    hdrdict={}

    for key in indices: 
        if key!='fdNUV' and key!='fdFUV': #Not full disc
            dat[key]=rasfits[indices[key]].data[...,lamwin[key][0]:lamwin[key][1]]*rcfs[key][None, None, :]

            dat[key][:,blanks[key]]=-200 #Reblanking the blanks

            hdrdict[key]=rasfits[indices[key]].header
            hdrdict[key]['CRVAL1']=wvlns[key][lamwin[key][0]]
            hdrdict[key]['NAXIS1']=lamwin[key][1]-lamwin[key][0]+1 #Counting "0", of course

            hdr0['TWMIN'+str(indices[key])]=wvlns[key][lamwin[key][0]]
            hdr0['TWMAX'+str(indices[key])]=wvlns[key][lamwin[key][1]]
            hdr0['TDMEAN'+str(indices[key])]=np.mean(dat[key])
            hdr0['TDRMS'+str(indices[key])]=np.sqrt(np.sum((dat[key]-np.mean(dat[key]))**2)/dat[key].size)
            hdr0['TDMEDN'+str(indices[key])]=np.median(dat[key])
            hdr0['TDMIN'+str(indices[key])]=np.min(dat[key])
            hdr0['TDMAX'+str(indices[key])]=np.max(dat[key])
            hdr0['TDSKEW'+str(indices[key])]=scist.skew(dat[key], axis=None)
            hdr0['TDKURT'+str(indices[key])]=scist.kurtosis(dat[key], axis=None)           
            p_results=np.percentile(dat[key], (1, 10, 25, 75, 90, 95, 98, 99))
            for entry, ptile in zip((1, 10, 25, 75, 90, 95, 98, 99), p_results):
                hdr0[f"TDP{entry:02d}_{str(indices[key])}"]=ptile
            
        else: #Full disc
            dat[key]=rasfits[indices[key]].data[lamwin[key][0]:lamwin[key][1]]*rcfs[key][:, None, None]

            hdrdict[key]=rasfits[indices[key]].header
            hdrdict[key]['CRVAL3']=wvlns[key][lamwin[key][0]]
            hdrdict[key]['NAXIS3']=lamwin[key][1]-lamwin[key][0]+1 #Counting "0", of course


    if key!='fdNUV' and key!='fdFUV': #Not full disc
        for ind, key in enumerate(dat):
            if ind==0:
                dattot=dat[key] #Needed for header stuff. (DATa TOTal)
            else:
                dattot=np.concatenate((dattot, dat[key]), axis=2)
        hdr0['DATAMEAN']=np.mean(dattot)
        hdr0['DATARMS']=np.sqrt(np.sum((dattot-np.mean(dattot))**2)/dattot.size)
        hdr0['DATAMEDN']=np.median(dattot)
        hdr0['DATAMIN']=np.min(dattot)
        hdr0['DATAMAX']=np.max(dattot)
        hdr0['DATASKEW']=scist.skew(dattot, axis=None)
        hdr0['DATAKURT']=scist.kurtosis(dattot, axis=None)

        p_results=np.percentile(dattot, (1, 10, 25, 75, 90, 95, 98, 99))
        for entry, ptile in zip((1, 10, 25, 75, 90, 95, 98, 99), p_results):
            hdr0[f"DATAP{entry:02d}"]=ptile

        del dattot 

        phdu=fits.PrimaryHDU(None, header=hdr0)
        hduls=[phdu]
        for key in indices:
            hduls.append(fits.ImageHDU(dat[key], header=hdrdict[key]))
    
    else:
        phdu=fits.PrimaryHDU(dat[list(indices.keys())[0]], header=hdrdict[list(indices.keys())[0]])
        hduls=[phdu]
        for hds in range(1, len(rasfits)):
            hduls.append(rasfits[hds])

    hdul=fits.HDUList(hduls)
    hdul.verify('fix')

    if error and ('fdNUV' not in indices):
        #Only for rasters currently
        #3.1 and 1.2 are the dark current values quoted in DN in "De Pontieu et al. 2014, Solar Physics, 289, 2733"
        dark={'FUV':3.1*response['DN2PHOT_SG'][FUVind], 'NUV':1.2*response['DN2PHOT_SG'][NUVind]}
        errs={}
        for i, key in enumerate(indices):
            temp=rasfits[indices[key]].data[...,lamwin[key][0]:lamwin[key][1]]
            if 'FUV' in rasfits[0].header['TDET'+str(indices[key])]:
                dd=dark['FUV']
                d2p=d2pFUV
            else:
                dd=dark['NUV']
                d2p=d2pNUV
            ph=temp*d2p 
            ph[ph<0]=0
            #Everything here is just a constant, so just multiply
            errs[i+1]=(np.sqrt(ph+dd**2))*np.abs(rcfs[key][None, None, :])

        hdr0e={i:rasfits[0].header[i] for i in rasfits[0].header['TDESC*']}
        phdue=fits.PrimaryHDU(None, header=fits.Header(hdr0e))
        hdulse=[phdue]
        hdrdicte=dc(hdrdict)
        for i in hdrdicte:
            for j in hdrdicte[i]['PC*']:
                del hdrdicte[i][j]
        for i, key in enumerate(indices):
            hdulse.append(fits.ImageHDU(errs[i+1], header=hdrdicte[key]))
        hdule=fits.HDUList(hdulse)
        hdule.verify('fix')  

    if error and ('fdNUV' not in indices):
        return(hdul, hdule)
    else:
        return(hdul)             
        