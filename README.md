# irispreppy
For radiometrically calibrating and PSF deconvolving IRIS data.

I dislike how I need to own proprietary software (IDL) just to simply prepare my data. I use Python for my analysis, why can't I radiometrically calibrate and deconvolve with it?
This has been a passion project of mine during my PhD. The radiometric calibration keeps itself up to date with the response files by checking https://hesperia.gsfc.nasa.gov/ssw/iris/response/ every time it is run. If it finds new files, it downloads them before continuing.

These scripts should be general purpose and "just work". No janky hacks are present.

This remains untested on Mac. However, I expect it to work on UNIX-like OSes.

---

Usage

irispreppy can take a list of HDU objects, a list of a directory paths to fits, or single HDU objects. For example,

```python
from astropy.io import fits
import irispreppy as ip

raw=fits.open("path/to/file.fits") #Raw data
rc=ip.radcal(raw)		   #Radiometrically calibrated
rc_d=ip.deconvolve(rc)		   #Radiometrically calibrated and deconvolved
```

---

Thank you to Dr Graham S. Kerr for IRIS_SG_deconvolve.py and IRIS_SG_PSFs.pkl.

Special thanks to Dr C.M.J. Osborne (https://github.com/Goobley) for putting up with my incessant and innane questions.
