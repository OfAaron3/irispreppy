import numpy as np
import scipy.fft as fft


def IRIS_SG_deconvolve(data_in, psf, 
                       iterations=10, 
                       fft_div=False,
                       ):


    '''

    Graham S. Kerr
    July 2020
 
    NAME:   IRIS_SG_Deconvolve.py

    PURPOSE: Deconvolves IRIS SG data using the PSFs from Courrier et al 2018.
    
    Input Parameters: 
        data_in    -- A 2D IRIS SG array [ypos, wavelength]
        psf        -- The appropriate PSF
        iterations -- The number of Richardson Lucy iterations to run through (Default = 10)
        fft_div    -- Set to use skip iterations and instead deconvolve by division Fourier Space   

    NOTES: Based on iris_sg_deconvolve.pro by Hans Courrier, but not all the functionality
           is included here yet

    
    History
    GSK 2020: Code translated
    '''


    #Remove negative values 
    dcvim = data_in.copy()
    dcvim[dcvim<0] = 0
    data_in_zr = dcvim

    if fft_div:
        dcvim = FFT_conv_1D(data_in,psf,div=True)
    else:
        for ind in range(1,iterations+1):
            step1 = data_in_zr/(FFT_conv_1D(dcvim,psf, rev_psf=False, div=False))
            step2 = FFT_conv_1D(step1,psf, rev_psf=True)
            dcvim = dcvim * step2

    return(dcvim)


def FFT_conv_1D(datain, psfin, div = False, rev_psf=False):
    
    '''   
    Notes AWP: This is not a good docstring, will fix later.

    Input parameters  
        datain -- a 2D data array [nominally, slit pos vs wavelength]
        psfin  -- the PSF to be applied in the y-direction
        div -- Set to True to divide in Fourier space (Default is False, so multiply in Fourier space)
        rev_psf -- Set to reverse the 1D input PSF
              
    Output:            
        dataout -- the input data convolved with the PSF
    
    Hisotory:
    GSK 2020: Pretty much copied exactly from Hans Courrier's IDL version in the SSW IRIS software tree, as part of iris_sg_deconvole.pro
    AWP 2021: Code fixed
    AWP 2023: Streamlined some of the code
    '''
  
    #length of input psf
    psflen = len(psfin)
    
    #dimensions of input data
    imsize = datain.shape
    
    #Get difference of image size and psf length
    ydiff = imsize[0]-psflen  
    
    #Cut the PSF if it is too long
    if ydiff < 0:
        pin=psfin[-(ydiff//2):ydiff//2+ydiff%2]
        #renormalize PSF (dx=1)
        pin = pin/np.sum(pin) 
        
    #Pad the PSF if it is too short
    if ydiff > 0:
        rs = ydiff//2
        pin=np.pad(psfin, [rs, rs+ydiff%2], 'edge')
        #Use edge values or you get a divide by zero error
        #renormalize PSF (dx=1)
        pin = pin/np.sum(pin) 
    
    #Replicate the PSF over wavelength array also
    pin_full = (np.repeat(pin, imsize[1])).reshape(len(pin), imsize[1])

    #Shift PSF center to zero to center output (slicing instead of rolling)
    pin_full = np.concatenate([pin_full[imsize[0]//2:], pin_full[:imsize[0]//2]])

    #Reverse the PSF if needed
    if rev_psf:
        pin_full = np.flip(pin_full,axis=0)        
        if psflen%2 == 0:
            pin_full = np.roll(pin_full,(1,0),axis=0)
    
    #Perform the FFT
    fpsf = fft.fft(pin_full,axis=0)
    datain=datain.astype(np.float64)
    fdatain = fft.fft(datain,axis=0)
    
    #Multiply(divide) the PSF and data, and convert back to k space
    if not div: 
        dataout = fft.ifft((fdatain*fpsf),axis=0).real
    else:
        dataout = fft.ifft((fdatain/fpsf),axis=0).real
    
    return(dataout)
