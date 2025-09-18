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
        Filename of output. If not set, will be saved as original filename (and path) with '_d' appended

    Example
    -------
    >>> from astropy.io import fits
    >>> import irispreppy as ip
    >>> f=fits.open('iris_raster.fits')
    >>> frc=ip.deconvolve_and_save(f)
    '''
    if filename is None:
        if '_rc'!=path.splitext(ras.filename())[0][-3:]:
            filename=path.splitext(ras.filename())[0]+'_d.fits'
        else:
            filename=path.splitext(ras.filename())[0]+'d.fits'
        
    hdul=deconvolve(ras)
    hudls.writeto(filename)
    