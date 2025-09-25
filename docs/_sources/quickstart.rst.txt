Quickstart
==========

Installation
------------

Installation should be simple through PyPI, and the package is pure python.

.. code-block:: bash

   python -m pip install irispreppy

Simplest Usage
--------------

Once installed, the simplest way to radiometrically calibrate and point spread function deconvolve your data is like so,

.. code-block:: python

    from astropy.io import fits
    import irispreppy as ip 

    f=fits.open("iris_raster.fits")
    frc=ip.radiometric_calibrate(f)
    frcd=ip.deconvolve(frc)

The files returned by irispreppy are structured identically to that of Level-2 IRIS fits. 

Saving files
------------
In order to radiometrically calibrate, point spread function deconvolve, and save your data you can do,

.. code-block:: python

    from astropy.io import fits
    import irispreppy as ip 

    f=fits.open("iris_raster.fits")
    ip.calibrate_and_save(f)
    frc=fits.open("iris_raster_rc.fits")
    ip.deconvolve_and_save(frc)
    frcd=fits.open("iris_raster_rcd.fits")

Where we also reopen the data that has been saved to disk.

`calibrate_and_save` will append `_rc` to the end of your filename, and `deconvolve` will append `_d` to your filename, unless it has `_rc` on the end already, where it will append `d`. To change this behaviour, simply specify the filename (and path) like so,

.. code-block:: python

    from astropy.io import fits
    import irispreppy as ip 

    f=fits.open("iris_raster.fits")
    ip.calibrate_and_save(f, filename="iris_raster_calibrated.fits")
    frc=fits.open("iris_raster_calibrated.fits")
    ip.deconvolve_and_save(frc, filename="iris_raster_calibrated_and_deconvolved.fits")
    frcd=fits.open("iris_raster_calibrated_and_deconvolved.fits"
