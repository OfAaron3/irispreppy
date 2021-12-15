from astropy.io import fits
import numpy as np
from glob import glob as ls
from IRIS_SG_deconvolve import *
import pickle

from copy import deepcopy as dc
import scipy.stats as scist
import concurrent.futures
from tqdm import tqdm

raspath='/home/aaron/Desktop/herc_dat/20191001/rc/*.fits'

rasdirec=ls(raspath)
rasdirec.sort()

with open('IRIS_SG_PSFs.pkl', 'rb') as psfpkl:
    psfsin=pickle.load(psfpkl)

psfs={'FUV1':psfsin['sg_psf_1336'], 'FUV2':psfsin['sg_psf_1394'], 'NUV':psfsin['sg_psf_2796']}


# for i in tqdm(range(0, len(rasdirec))):
#     rass=rasdirec[i]
def ParDecon(rass, psfs):
    rasfits=fits.open(rass)
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
            decondict[key][j]=IRIS_SG_deconvolve(rasfits[indices[key]].data[j], psf=psfs[psfind], fft_div=True)


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

    if '_rc' in rass:
        hdul.writeto(rass[:-8]+'_rd.fits')
    else:
        hdul.writeto(rass[:-5]+'_decon.fits')


with concurrent.futures.ProcessPoolExecutor() as executor:
    futures=[executor.submit(ParDecon, rass=rasdirec[i], psfs=psfs) for i in range(0, len(rasdirec))]
    for f in tqdm(concurrent.futures.as_completed(futures), total=len(rasdirec)):
        pass
for f in enumerate(futures):
    test=f[1].result() #Just to make sure nothing crashed