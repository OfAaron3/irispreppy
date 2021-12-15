# IRISPreppy
For radiometrically calibrating and PSF deconvolving IRIS data.

I dislike how I need to own proprietary software (IDL) just to simply prepare my data. I use Python for my analysis, why can't I radiometrically calibrate and deconvolve with it?
This has been a passion project of mine during my PhD. It keeps itself up to date with the response files by checking https://hesperia.gsfc.nasa.gov/ssw/iris/response/ every time it is run. If it finds new files, it downloads them before continuing.

These scripts should be general purpose and "just work". No janky hacks are present.

Just need to write some scripts as functions and we're good to go.

~~~

IRIS_SG_deconvolve.py authored by Dr Graham S. Kerr.
IRIS_SG_PSFs.pkl supplied by Dr Graham S. Kerr.

iris_get_response.pro translated by Aaron W. Peat.
fit_iris_xput.pro translated by Aaron W. Peat.
radcal.py authored by Aaron W. Peat
deconvolve.py authored by Aaron W. Peat

~~~

Special thanks to Dr C.M. Osborne (https://github.com/Goobley) for putting up with my incessant and innane questions.
