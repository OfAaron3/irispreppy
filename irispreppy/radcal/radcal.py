import datetime as dt
from copy import deepcopy as dc
from glob import glob as ls
from subprocess import call as bashrun

import numpy as np
import scipy.stats as scist
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.io import readsav
from tqdm import tqdm

from . import iris_get_response as igr

#################################################################

def radcal(ras, save=False, quiet=True, debug=False):
    '''Radiometric Calibration of IRIS files
    Input Paramaters:
        ras: String, list, or astropy.io.fits.hdu.hdulist.HDUList (hdu)
             String: Path to IRIS spectrograph file
                     Path to IRIS files using wildcard (ie, /path/to/files/*fits) of same observation
             List: List of paths to spectrograph file of same observation
                   List of hdus from same observation
             hdu : An IRIS observation
        save: Save the output files in the same directory as the input file with "rc" appended to the end of the file. Default is False
              If True: Saves radiometrically calibrated files in the same directory as the input file with "rc" appended to the end of the filename. Does not return radiometrically calibrated hduls
              If False: Does exact opposite of True case. Does not save, and does return radiometrically calibrated hduls
        quiet: Suppress print outputs. Default True.
    '''

    pathlistin=False
    hdulistin=False
    if type(ras)==fits.hdu.hdulist.HDUList:
        assert ras[0].header['TELESCOP']=='IRIS'
        rasfits=dc(ras)

    elif '*' in ras:
        ras=ls(rass)
        ras.sort()
        assert fits.open(ras[0]).header['TELESCOP']=='IRIS'
        pathlistin=True

    elif type(ras)==str:
        try:
            rasfits=fits.open(ras)
            assert rasfits[0].header['TELESCOP']=='IRIS'
        except NameError:
            raise ValueError("Must supply fits file or path to fits file or * directory for one set of observations")

    elif type(ras)==list:
        if type(ras[0])==fits.hdu.hdulist.HDUList:
            assert ras[0].header['TELESCOP']=='IRIS'
            hdulistin=True
        else:
            try:
                assert fits.open(ras[0])[0].header['TELESCOP']=='IRIS'
                pathlistin=True
            except NameError:
                raise ValueError("Must supply fits file or * directory for one set of observations")
    else:
        raise ValueError("Must supply fits file or * directory for one set of observations")


    ############
    # Response #
    ############
    if pathlistin:
        rasfits=fits.open(ras[0])
    elif hdulistin:
        rasfits=ras[0]

    begin=dt.datetime.strptime(rasfits[0].header['STARTOBS'], '%Y-%m-%dT%H:%M:%S.%f')
    end=dt.datetime.strptime(rasfits[0].header['ENDOBS'], '%Y-%m-%dT%H:%M:%S.%f')
    midtime=dt.datetime.strftime((begin+((end-begin)/2)), '%Y-%m-%dT%H:%M:%S.%fZ')

    response=(igr.iris_get_response(midtime))[0]

    ###################################
    # Find which index is FUV and NUV #
    ###################################
    FUVind=-1
    NUVind=-1
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

    ###################
    # Wavelength axes #
    ###################

    FUV=response['LAMBDA']*(response['AREA_SG'][FUVind]>0).astype(int)
    flag=False
    FUVcutoff=-1

    while FUV[FUVcutoff]!=0 or not flag:
        FUVcutoff+=1
        if FUV[FUVcutoff]>0 and not flag:
            flag=True

    FUV1=FUV[:FUVcutoff]
    FUV1=(FUV1[FUV1!=0])*10

    FUV2=FUV[FUVcutoff:]
    FUV2=(FUV2[FUV2!=0])*10

    del flag, FUV

    NUV=response['LAMBDA']*(response['AREA_SG'][NUVind]>0).astype(int)
    NUV=(NUV[NUV!=0])*10

    #################
    # Photon Energy #
    #################
    #Lambda is in nm. 1e7 ergs = 1 Joule
    h=6.62607004e-34
    c=3e8
    eFUV1=1e7*h*c/(FUV1*10e-10)
    eFUV2=1e7*h*c/(FUV2*10e-10)
    eNUV=1e7*h*c/(NUV*10e-10)

    ##################
    # Effective Area #
    ##################
    aFUV1=(response['AREA_SG'][FUVind][:FUVcutoff])[response['AREA_SG'][FUVind][:FUVcutoff]>0]
    aFUV2=(response['AREA_SG'][FUVind][FUVcutoff:])[response['AREA_SG'][FUVind][FUVcutoff:]>0]
    aNUV=response['AREA_SG'][NUVind][response['AREA_SG'][NUVind]>0]

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
        if rasfits[0].header['TDET'+str(indices[key])]=='FUV':
            const[key]=d2pFUV/(pixxy[key]*pixl[key]*tfuv*wslit)
        else:
            const[key]=d2pNUV/(pixxy[key]*pixl[key]*tnuv*wslit)

    ###############################################
    # Interpolating wavelength dependant function #
    ###############################################

    #f(\lambda)=E(\lambda)/Aeff(\lambda)
    f_FUV1=interp1d(FUV1, eFUV1/aFUV1, kind='cubic')
    f_FUV2=interp1d(FUV2, eFUV2/aFUV2, kind='cubic')
    f_NUV=interp1d(NUV, eNUV/aNUV, kind='cubic')


    ###########################################################
    # Wavelength Trimming and Radiometric Calibration Factors #
    ###########################################################
    
    lamwin={} #LAMbda WINdow
    wvlns={}
    rcfs={} #Radiometric Calibration FactorS
    for key in indices:
        if rasfits[0].header['TDET'+str(indices[key])]=='FUV1':
            wvlns[key]=np.add(np.multiply(np.subtract(np.arange(0, rasfits[indices[key]].header['NAXIS1']), rasfits[indices[key]].header['CRPIX1']-1), rasfits[indices[key]].header['CDELT1']), rasfits[indices[key]].header['CRVAL1'])
            lamwin[key]=[-1, -1]
            for ind, wvln in enumerate(wvlns[key]):
                if wvln>=FUV1[0] and lamwin[key][0]==-1:
                    lamwin[key][0]=ind
                elif wvln>FUV1[-1] and lamwin[key][1]==-1:
                    lamwin[key][1]=ind-1
                    break
            rcfs[key]=f_FUV1(wvlns[key][lamwin[key][0]:lamwin[key][1]])*const[key]

        elif rasfits[0].header['TDET'+str(indices[key])]=='FUV2':
            wvlns[key]=np.add(np.multiply(np.subtract(np.arange(0, rasfits[indices[key]].header['NAXIS1']), rasfits[indices[key]].header['CRPIX1']-1), rasfits[indices[key]].header['CDELT1']), rasfits[indices[key]].header['CRVAL1'])
            lamwin[key]=[-1, -1]
            for ind, wvln in enumerate(wvlns[key]):
                if wvln>=FUV2[0] and lamwin[key][0]==-1:
                    lamwin[key][0]=ind
                elif wvln>FUV2[-1] and lamwin[key][1]==-1:
                    lamwin[key][1]=ind-1
                    break
            rcfs[key]=f_FUV2(wvlns[key][lamwin[key][0]:lamwin[key][1]])*const[key]

        elif rasfits[0].header['TDET'+str(indices[key])]=='NUV':
            wvlns[key]=np.add(np.multiply(np.subtract(np.arange(0, rasfits[indices[key]].header['NAXIS1']), rasfits[indices[key]].header['CRPIX1']-1), rasfits[indices[key]].header['CDELT1']), rasfits[indices[key]].header['CRVAL1'])
            lamwin[key]=[-1, -1]
            for ind, wvln in enumerate(wvlns[key]):
                if wvln>=NUV[0] and lamwin[key][0]==-1:
                    lamwin[key][0]=ind
                elif wvln>NUV[-1] and lamwin[key][1]==-1:
                    lamwin[key][1]=ind-1
                    break
            rcfs[key]=f_NUV(wvlns[key][lamwin[key][0]:lamwin[key][1]])*const[key]

        else:
            raise Exception("You have detectors that are not FUV1, FUV2, or NUV in your fits file.")

    ###########################
    # Radiometric Calibration #
    ###########################
    if not quiet:
        print("Creating and saving calibrated spectrograph fits...")    
        print("(This may take some time as it needs to recalculate the excessive amount of stats in the headers)")

    if hdulistin or pathlistin:
        end=len(ras)
    else:
        end=1
    
    for k in tqdm(range(0, end), disable=quiet):
        if pathlistin:
            rasfits=fits.open(ras[k])
        elif hdulistin:
            rasfits=ras[k]

        hdr0=dc(rasfits[0].header)
        hdr0['HISTORY']='Mgii, SiIV and CII radiometric calibration performed on '+dt.datetime.now().strftime("%Y-%m-%d")
        hdr0['HISTORY']='FITS made with astropy on '+dt.datetime.now().strftime("%Y-%m-%d")
        hdr0['BUNIT']="erg s^-1 cm^-2 angstrom^-1 sr^-1"
        dat={}
        hdrdict={}

        for key in indices: 
            dat[key]=dc(rasfits[indices[key]].data[...,lamwin[key][0]:lamwin[key][1]])
            for ind, _ in np.ndenumerate(dat[key][...,0]):
                dat[key][ind]=np.multiply(dat[key][ind], rcfs[key])
            hdrdict[key]=dc(rasfits[indices[key]].header)
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
            
            flatdat=np.sort(dat[key].flatten())
            hdr0['TDP01_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.01))]
            hdr0['TDP10_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.1))]
            hdr0['TDP25_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.25))]
            hdr0['TDP75_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.75))]
            hdr0['TDP90_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.9))]
            hdr0['TDP95_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.95))]
            hdr0['TDP98_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.98))]
            hdr0['TDP99_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.99))]
            del flatdat

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
        flatdattot=np.sort(dattot.flatten())
        hdr0['DATAP01']=flatdattot[int(np.round(len(flatdattot)*0.01))]
        hdr0['DATAP10']=flatdattot[int(np.round(len(flatdattot)*0.1))]
        hdr0['DATAP25']=flatdattot[int(np.round(len(flatdattot)*0.25))]
        hdr0['DATAP75']=flatdattot[int(np.round(len(flatdattot)*0.75))]
        hdr0['DATAP90']=flatdattot[int(np.round(len(flatdattot)*0.9))]
        hdr0['DATAP95']=flatdattot[int(np.round(len(flatdattot)*0.95))]
        hdr0['DATAP98']=flatdattot[int(np.round(len(flatdattot)*0.98))]
        hdr0['DATAP99']=flatdattot[int(np.round(len(flatdattot)*0.99))]
        del dattot, flatdattot #I imagine these are large, so delete them after they are no longer needed

        phdu=fits.PrimaryHDU(None, header=hdr0)
        hduls=[phdu]
        for key in indices:
            hduls.append(fits.ImageHDU(dat[key], header=hdrdict[key]))
        hdul=fits.HDUList(hduls)
        if save:
            hdul.writeto(rasfits.filename()[:-5]+'_rc.fits')
        else:
            if hdulistin or pathlistin:
                if k!=0:
                    callist.append(hdul)
                elif k==0:
                    callist=[hdul]                
            
    if not quiet:
        print("Done!")

    if hdulistin or pathlistin:
        return(callist)
    else:
        return(hdul)
