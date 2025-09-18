.. irispreppy documentation master file, created by
   sphinx-quickstart on Wed Sep 17 16:48:21 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to irispreppy's documentation
=====================================


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   irispreppyapi
   

Introduction
------------

irispreppy is built to perform radiometric calibration and point spread function deconvolution of data from the spectrograph onboard the Interface Region Imaging Spectrograph (IRIS_). 

The radiometric calibration follows the procedure laid out in section 5.2 of IRIS Technical Notes 26 (ITN 26; `html <https://iris.lmsal.com/itn26/calibration.html#radiometric-calibration>`_; `pdf <http://iris.lmsal.com/itn26/itn26.pdf>`_), and more generally IRIS Technical Notes 24 (ITN 24; `pdf <https://www.lmsal.com/iris_science/doc?cmd=dcur&proj_num=IS0123&file_type=pdf>`_).

The point spread function deconvolution is done using either the Richardson_-Lucy_ deconvolution method, or through division in Fourier space. The number of default Richardson_-Lucy_ iterations used is discussed in `Courrier et al. (2018) <https://doig.org/10.1007/s11207-018-1347-9>`_.

.. _IRIS: https://doi.org/10.1007/s11207-014-0485-y 
.. _Richardson: https://doi.org/10.1364/JOSA.62.000055
.. _Lucy: https://doi.org/10.1086/111605



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

