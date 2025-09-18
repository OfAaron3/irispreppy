from .radcal import radiometric_calibrate as radcal
from os import path

def calibrate_and_save(ras, errors=False, filename=None):
    '''
    Performs radiometric calibration of IRIS spectrograph files using default values and saves.\n
    (See Section 5.2 of ITN26 for more details - https://iris.lmsal.com/itn26/calibration.html#radiometric-calibration)
    
    Parameters
    ----------
    ras : astropy.io.fits.hdu.hdulist.HDUList
        Input IRIS raster
    error : bool
        Whether to calculate errors (beta). Ignored if `ras` is a full disc mosaic. Default: False
    filename : string
        Filename of output. If not set, will be saved as original filename (and path) with '_rc' appended

    Example
    -------
    >>> from astropy.io import fits
    >>> import irispreppy as ip
    >>> f=fits.open('iris_raster.fits')
    >>> frc=ip.calibrate_and_save(f)
    '''
    if filename is None:
        filename=path.splitext(ras.filename())[0]+'_rc.fits'
    if errors:
        filenamee=path.splitext(filename)[0]+'e.fits'
        
    hduls=radcal(ras, error=error)
    if type(hduls)==tuple: 
        #Error is ignored for full disc mosaics, have to check what is returned
        hdul=hduls[0]
        hdule=hduls[1]
        hdul.writeto(filename)
        hdule.writeto(filenamee)
    else:
        if error:
            import warnings
            warnings.warn('Errors are currently not available for full disc mosaics.')
        hudls.writeto(filename)
