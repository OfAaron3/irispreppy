import concurrent.futures
import pickle
from copy import deepcopy as dc
from glob import glob as ls
from os import cpu_count as cpus
from os import path

import numpy as np
import scipy.stats as scist
from astropy.io import fits
from tqdm import tqdm

from . import IRIS_SG_deconvolve as isd

#There are two functions here


def ParDecon(rasfits, psfs, save=False):
    '''Function acts as a wrapper around IRIS_SG_deconvolve
    Input Paramteres: 
        rasfits: The hdu to be deconvolved
        psfs: The point spread functions in a dictionary
        save: If True: Save the files with d appended
              If False: Return the deconvolved hdus
    Output:
        If save=False: Deconcolved hdu
        If save=True: 0
    '''
    hdr0=dc(rasfits[0].header)
    nlines=rasfits[0].header['NEXP']
    indices={rasfits[0].header[name]: ind+1 for ind, name in enumerate(rasfits[0].header['TDESC*'])}
    decondict={}
    hdrdict={}
    for key in indices:
        decondict[key]=np.zeros_like(rasfits[indices[key]].data)
        psfind=rasfits[0].header['TDET'+str(indices[key])]
        hdrdict[key]=dc(rasfits[indices[key]].header)
        for j in range(0, nlines):
            decondict[key][j]=isd.IRIS_SG_deconvolve(rasfits[indices[key]].data[j], psf=psfs[psfind], fft_div=True)


        hdr0['TDMEAN'+str(indices[key])]=np.mean(decondict[key])
        hdr0['TDRMS'+str(indices[key])]=np.sqrt(np.sum((decondict[key]-np.mean(decondict[key]))**2)/decondict[key].size)
        hdr0['TDMEDN'+str(indices[key])]=np.median(decondict[key])
        hdr0['TDMIN'+str(indices[key])]=np.min(decondict[key])
        hdr0['TDMAX'+str(indices[key])]=np.max(decondict[key])
        hdr0['TDSKEW'+str(indices[key])]=scist.skew(decondict[key], axis=None)
        hdr0['TDKURT'+str(indices[key])]=scist.kurtosis(decondict[key], axis=None)

        flatdat=np.sort(decondict[key].flatten())
        hdr0['TDP01_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.01))]
        hdr0['TDP10_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.1))]
        hdr0['TDP25_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.25))]
        hdr0['TDP75_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.75))]
        hdr0['TDP90_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.9))]
        hdr0['TDP95_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.95))]
        hdr0['TDP98_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.98))]
        hdr0['TDP99_'+str(indices[key])]=flatdat[int(np.round(len(flatdat)*0.99))]
        del flatdat

    for ind, key in enumerate(decondict):
        if ind==0:
            dattot=decondict[key] #Needed for header stuff. (DATa TOTal)
        else:
            dattot=np.concatenate((dattot, decondict[key]), axis=2)
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
        hduls.append(fits.ImageHDU(decondict[key], header=hdrdict[key]))
    hdul=fits.HDUList(hduls)  

    if save:
        hdul.writeto(rasfits.filename()[:-5]+'d.fits')
        return(0)
    else:
        return(hdul)


def deconvolve(ras, quiet=False, save=False, limitcores=False):
    '''Function prepares input to ParDecon
    Input Paramteres: 
        ras: String, list, or astropy.io.fits.hdu.hdulist.HDUList (hdu)
             String: Path to IRIS spectrograph file
                     Path to IRIS files using wildcard (ie, /path/to/files/*fits) of same observation
             List: List of paths to spectrograph file of same observation
                   List of hdus from same observation
             hdu : An IRIS observation
        quiet: If True, suppress all print statements
        save: If True: Save the files with d appended
              If False: Return the deconvolved hdus
        limitcores: If True: use all but one core. If False use all cores. 
    Output:
        If save=False: Deconcolved hdu(s). 
        If save=True: 0
    '''

    nworkers=cpus()-int(limitcores)
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

    toppath=path.dirname(path.realpath(__file__))
    with open(toppath+'/IRIS_SG_PSFs.pkl', 'rb') as psfpkl:
        psfsin=pickle.load(psfpkl)
    
    psfs={'FUV1':psfsin['sg_psf_1336'], 'FUV2':psfsin['sg_psf_1394'], 'NUV':psfsin['sg_psf_2796']}      

    if pathlistin:
        with concurrent.futures.ProcessPoolExecutor(workers=nworkers) as executor:
            futures=[executor.submit(ParDecon, rasfits=fits.open(ras[i]), psfs=psfs, save=save) for i in range(0, len(ras))]
            for f in tqdm(concurrent.futures.as_completed(futures), total=len(rasdirec), disable=quiet):
                pass 
        out=[f for f in futures]
    elif hdulistin:
        with concurrent.futures.ProcessPoolExecutor(workers=nworkers) as executor:
            futures=[executor.submit(ParDecon, rasfits=ras[i], psfs=psfs, save=save) for i in range(0, len(ras))]
            for f in tqdm(concurrent.futures.as_completed(futures), total=len(rasdirec), disable=quiet):
                pass
        out=[f for f in futures]
    else:
        out=ParDecon(rasfits=ras, psfs=psfs, save=save)
    
    if not save:
        return(out)
    else:
        return(0)
