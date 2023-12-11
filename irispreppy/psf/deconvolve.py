import pickle
from copy import deepcopy as dc
from os import path

import numpy as np
import scipy.stats as scist
from astropy.io import fits

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
    nlines=hdr0['NEXP']
    indices={hdr0[name]: ind+1 for ind, name in enumerate(hdr0['TDESC*'])}
    deconlst=[]
    for index, key in enumerate(indices):
        deconlst.append(np.zeros_like(rasfits[indices[key]].data))
        psfind=hdr0['TDET'+str(indices[key])]
        for j in range(0, nlines):
            deconlst[index][j]=isd.IRIS_SG_deconvolve(rasfits[indices[key]].data[j], psf=psfs[psfind], fft_div=True)

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

    phdu=fits.PrimaryHDU(None, header=hdr0)
    hduls=[phdu]
    for index, key in enumerate(indices):
        hduls.append(fits.ImageHDU(deconlst[index], header=rasfits[indices[key]].header))
    hdul=fits.HDUList(hduls)  

    if save:
        hdul.writeto(path.splitext(rasfits.filename())[0]+'d.fits')
        return(0)
    else:
        return(hdul)


def deconvolve(ras, save=False):
    '''Function prepares input to ParDecon
    Input Paramteres: 
        ras: astropy.io.fits.hdu.hdulist.HDUList of an IRIS observation
        save: If True: Save the files with d appended
              If False: Return the deconvolved hdus
        limitcores: If True: use all but one core. If False use all cores. 
    Output:
        If save=False: Deconcolved hdu(s). 
        If save=True: 0
    '''

    toppath=path.dirname(path.realpath(__file__))
    with open(toppath+'/IRIS_SG_PSFs.pkl', 'rb') as psfpkl:
        psfsin=pickle.load(psfpkl)
    
    psfs={'FUV1':psfsin['sg_psf_1336'], 'FUV2':psfsin['sg_psf_1394'], 'NUV':psfsin['sg_psf_2796']}      

    out=ParDecon(rasfits=ras, psfs=psfs, save=save)
    
    if not save:
        return(out)
    else:
        return(0)
