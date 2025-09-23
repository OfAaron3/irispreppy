from .deconvolve import deconvolve
from os import path

def deconvolve_and_save(ras, filename=None):
    '''
    Deconvolves IRIS PSF from IRIS spectrograph files using default values and saves.\n
    (See Courrier et al. 2018 for more information - DOI: 10.1007/s11207-018-1347-9)
    
    Parameters
    ----------
    ras : astropy.io.fits.hdu.hdulist.HDUList
        Input IRIS raster

    filename : string
        Filename of output. If not set, will be saved as original filename (and path) with '_d' appended.

    Notes
    -----
    If no filename is present in the HDUL and filename is not set, the original filename will be deciphered from the header information.

    Example
    -------
    >>> from astropy.io import fits
    >>> import irispreppy as ip
    >>> f=fits.open('iris_raster.fits')
    >>> frc=ip.deconvolve_and_save(f)
    '''
    if filename is None:
        if ras.filename()!=None:
            if path.splitext(ras.filename())[0][-3:]!='_rc':
                filename=path.splitext(ras.filename())[0]+'_d.fits'
            else:
                filename=path.splitext(ras.filename())[0]+'d.fits'
        else:
            hdr=ras[0].header
            filename='./iris_l'+str(int(hdr['DATA_LEV']))+'_'+hdr['DATE_OBS'].replace('-','').replace(':', '').replace('T', '_')[:-4]+'_'+hdr['OBSID']+'_raster_t000_r'+str(hdr['RASRPT']-1).zfill(5)+'_d.fits'
            print("Filename not present in HDUL and no filename set. Saving to,")
            print(filename)

    hdul=deconvolve(ras)
    hdul.writeto(filename)
    