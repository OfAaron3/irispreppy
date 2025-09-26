import datetime as dt
import pickle
from copy import deepcopy as dc
from os import path

import numpy as np
import scipy.stats as scist
from astropy.io import fits

from . import IRIS_SG_deconvolve as isd


def decon(rasfits, psfs, iterations=0, fft=False):
    '''
    Wrapper around IRIS_SG_deconvolve

    Parameters
    ----------
        rasfits : astropy.io.fits.hdu.hdulist.HDUList
            The hdu to be deconvolved
        psfs : dict
            The point spread functions
        iterations : int or list
            How many iterations to use in Richardson-Lucy deconvolution. Integer input will be applied to all ImageHDUs. List input should be a separate number of iterations for each ImageHDU. Defaults - FUV: 10, NUV: 50
        fft : bool
            Whether to deconvolve by division in Fourier space instead of using a Richardson-Lucy deconvolution. Default: False
            
    Returns
    -------
        hdul : astropy.io.fits.hdu.hdulist.HDUList
            Deconvolved hdu
    '''
    hdr0=dc(rasfits[0].header)
    nlines=hdr0['NEXP']
    indices={hdr0[name]: ind+1 for ind, name in enumerate(hdr0['TDESC*'])}
    deconlst=[]
    for index, key in enumerate(indices):
        deconlst.append(np.zeros_like(rasfits[indices[key]].data))
        psfind=hdr0['TDET'+str(indices[key])]
        sbf=int(np.round(rasfits[indices[key]].header['CDELT2']*6))
        if type(iterations)==int:
            if iterations==0 and not fft:
                if psfind=='NUV':
                    its=50
                else:
                    its=10
            else:
                its=iterations
        elif type(iterations)==list or type(iterations)==np.ndarray:
            if len(iterations)!=len(indices):
                ValueError('List or arrays of iterations must be same length as number of wavelength windows and in same order as in fits file.')
            its=iterations[index]
        if type(its)!=int:
            if int(its)==its:
                its=int(its)
            else:
                raise TypeError("Number of iterations must be an integer (or list of integers).")
        if sbf!=1:
            psf=psfs[psfind]
            psf2use=np.array([np.sum(psf[i*sbf:i*sbf+sbf]) for i in range(0, int(len(psf)/sbf))])
        else:
            psf2use=psfs[psfind]
        for j in range(0, nlines):
            deconlst[index][j]=isd.IRIS_SG_deconvolve(rasfits[indices[key]].data[j], psf=psf2use, iterations=its, fft_div=fft)
        if not fft:
            deconlst[index][np.isnan(deconlst[index])]=0

        hdr0['TDMEAN'+str(indices[key])]=np.mean(deconlst[index])
        hdr0['TDRMS'+str(indices[key])]=np.sqrt(np.sum((deconlst[index]-np.mean(deconlst[index]))**2)/deconlst[index].size)
        hdr0['TDMEDN'+str(indices[key])]=np.median(deconlst[index])
        hdr0['TDMIN'+str(indices[key])]=np.min(deconlst[index])
        hdr0['TDMAX'+str(indices[key])]=np.max(deconlst[index])
        hdr0['TDSKEW'+str(indices[key])]=scist.skew(deconlst[index], axis=None)
        hdr0['TDKURT'+str(indices[key])]=scist.kurtosis(deconlst[index], axis=None)

        hdr0['TDP01_'+str(indices[key])]=np.percentile(deconlst[index], 1)
        hdr0['TDP10_'+str(indices[key])]=np.percentile(deconlst[index], 10)
        hdr0['TDP25_'+str(indices[key])]=np.percentile(deconlst[index], 25)
        hdr0['TDP75_'+str(indices[key])]=np.percentile(deconlst[index], 75)
        hdr0['TDP90_'+str(indices[key])]=np.percentile(deconlst[index], 90)
        hdr0['TDP95_'+str(indices[key])]=np.percentile(deconlst[index], 95)
        hdr0['TDP98_'+str(indices[key])]=np.percentile(deconlst[index], 98)
        hdr0['TDP99_'+str(indices[key])]=np.percentile(deconlst[index], 99)

    for ind, _ in enumerate(deconlst):
        if ind==0:
            dattot=deconlst[ind] #Needed for header stuff. (DATa TOTal)
        else:
            dattot=np.concatenate((dattot, deconlst[ind]), axis=2)
    hdr0['DATAMEAN']=np.mean(dattot)
    hdr0['DATARMS']=np.sqrt(np.sum((dattot-np.mean(dattot))**2)/dattot.size)
    hdr0['DATAMEDN']=np.median(dattot)
    hdr0['DATAMIN']=np.min(dattot)
    hdr0['DATAMAX']=np.max(dattot)
    hdr0['DATASKEW']=scist.skew(dattot, axis=None)
    hdr0['DATAKURT']=scist.kurtosis(dattot, axis=None)

    hdr0['DATAP01']=np.percentile(dattot, 1)
    hdr0['DATAP10']=np.percentile(dattot, 10)
    hdr0['DATAP25']=np.percentile(dattot, 25)
    hdr0['DATAP75']=np.percentile(dattot, 75)
    hdr0['DATAP90']=np.percentile(dattot, 90)
    hdr0['DATAP95']=np.percentile(dattot, 95)
    hdr0['DATAP98']=np.percentile(dattot, 98)
    hdr0['DATAP99']=np.percentile(dattot, 99)
    del dattot  

    if fft:
        hdr0['HISTORY']='PSF deconvolution performed by division in Fourier space on '+dt.datetime.now().strftime("%Y-%m-%d")
    else:
        hdr0['HISTORY']='PSF deconvolution performed by Richardson-Lucy algorithm on '+dt.datetime.now().strftime("%Y-%m-%d")
    hdr0['HISTORY']='FITS made with astropy on '+dt.datetime.now().strftime("%Y-%m-%d")

    phdu=fits.PrimaryHDU(None, header=hdr0)
    hduls=[phdu]
    for index, key in enumerate(indices):
        for yy in range(0, rasfits[indices[key]].header['NAXIS2']):
            if (rasfits[indices[key]].data[:,yy]<0).all():
                deconlst[index][:,yy]=-200 #Reblanking the blanks
        hduls.append(fits.ImageHDU(deconlst[index], header=rasfits[indices[key]].header))
    hdul=fits.HDUList(hduls)  

    return(hdul)


def deconvolve(ras, iterations=0, fft=False):
    '''
    Deconvolves IRIS PSF from IRIS spectrograph files.\n
    (See Courrier et al. 2018 for more information - DOI: 10.1007/s11207-018-1347-9)
    
    Parameters
    ----------
    ras : astropy.io.fits.hdu.hdulist.HDUList
        Input IRIS raster
    iterations : int or list
        How many iterations to use in Richardson-Lucy deconvolution. Integer input will be applied to all ImageHDUs. List input should be a separate number of iterations for each ImageHDU. Defaults - FUV: 10, NUV: 50
    fft : bool
        Whether to deconvolve by division in Fourier space instead of using a Richardson-Lucy deconvolution. Default: False

    Returns
    -------
    hdul : astropy.io.fits.hdu.hdulist.HDUList
        The deconvolved data

    Example
    -------
    >>> from astropy.io import fits
    >>> import irispreppy as ip
    >>> f=fits.open('iris_raster.fits')
    >>> frc=ip.deconvolve(f)

    '''

    if ras[0].header['NAXIS']!=0: #FD mosaic
        raise ValueError("PSF deconvolution of full disc mosaics is not possible.")
    toppath=path.dirname(path.realpath(__file__))
    with open(toppath+'/IRIS_SG_PSFs.pkl', 'rb') as psfpkl:
        psfsin=pickle.load(psfpkl)
    
    if ras[0].header['HISTORY'][-1][:-11]=='FITS made with astropy on':
        from datetime import datetime as dt
        import warnings
        if dt.strptime(ras[0].header['HISTORY'][-1][-10:], '%Y-%m-%d')<dt.strptime('2025-04-29', '%Y-%m-%d'):
            warnings.warn("Future versions of deconvolve.py may behave unexpectedly with rasters calibrated before 2025-04-29. If you plan to deconvolve again, please calibrate again", UserWarning)
            

    psfs={'FUV1':psfsin['sg_psf_1336'], 'FUV2':psfsin['sg_psf_1394'], 'NUV':psfsin['sg_psf_2796']}      

    out=decon(rasfits=ras, psfs=psfs, iterations=iterations, fft=fft)
    
    return(out)
