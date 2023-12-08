# irispreppy
For radiometrically calibrating and PSF deconvolving IRIS data.

To install run `pip install irispreppy`.

I dislike how I need to own proprietary software (IDL) just to simply prepare my data. I use Python for my analysis, why can't I radiometrically calibrate and deconvolve with it?
This has been a passion project of mine during my PhD. The radiometric calibration keeps itself up to date with the response files by checking https://hesperia.gsfc.nasa.gov/ssw/iris/response/ every time it is run. If it finds new files, it downloads them before continuing.

These scripts should be general purpose and "just work". No janky hacks are present.

This remains untested on Mac. However, I expect it to work on UNIX-like OSes.

---

tl;dr usage

irispreppy takes a single HDU objects. To calibrate and deconvolve,

```python
from astropy.io import fits
import irispreppy as ip

raw=fits.open("path/to/file.fits") #Raw data
rc=ip.radcal(raw)		   #Radiometrically calibrated
rc_d=ip.deconvolve(rc)		   #Radiometrically calibrated and deconvolved
```

To calibrate and deconvolve, and save,

```python
from astropy.io import fits
import irispreppy as ip

raw=fits.open("path/to/file.fits") #Raw data
ip.radcal(raw, save=True)		   #Radiometrically calibrated
rc=fits.open("path/to/file_rc.fits")
ip.deconvolve(rc, save=True)		   #Radiometrically calibrated and deconvolved
```

---

More in depth documentation will be added in the future.

---

Acknowledgements

Thank you to [Dr Graham S. Kerr](https://github.com/grahamkerr) for IRIS_SG_deconvolve.py and IRIS_SG_PSFs.pkl.

Special thanks to [Dr C.M.J. Osborne](https://github.com/Goobley) for putting up with my incessant and innane questions.

Makes use of the excellent WENO4 algorithm ([Janett et al. 2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.104J/abstract)) implemented in Python3 [by Dr C.M.J. Osborne here](https://github.com/Goobley/Weno4Interpolation).
