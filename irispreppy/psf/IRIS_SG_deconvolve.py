import numpy as np
import scipy.fft as fft
import warnings


def IRIS_SG_deconvolve(data_in, psf, iterations, fft_div=False):
    '''
    Deconvolves IRIS point spread function from IRIS spectrograph slits 
    (See Courrier et al. 2018 for more information - DOI: 10.1007/s11207-018-1347-9)
    
    Parameters
    ----------
    data_in : array_like
        A 2D IRIS SG array (a single slit position)
    psf : array_like
        The appropriate point spread function for the considered detector
    iterations : int 
        The number of Richardson Lucy iterations 
    fft_div : bool 
        Whether to deconvolve by division in Fourier space instead of using a Richardson-Lucy deconvolution. Default: False  

    Returns
    -------
    dcvim : array_like
        Deconvolved IRIS slit
    '''

    #Remove negative values 
    dcvim = data_in.copy()
    dcvim[dcvim<0] = 0
    data_in_zr = dcvim

    if fft_div:
        dcvim = FFT_conv_1D(data_in, psf, div=True)
    else:
        with warnings.catch_warnings(): #Divide by zero warnings will happen
            warnings.simplefilter('ignore')
            for ind in range(0,iterations):
                step1=data_in_zr/(FFT_conv_1D(dcvim, psf, rev_psf=False, div=False))
                step2=FFT_conv_1D(step1, psf, rev_psf=True)
                dcvim=dcvim*step2

    return(dcvim)


def FFT_conv_1D(datain, psfin, div=False, rev_psf=False):
    
    '''   
    Fourier Transform 1D Convolution

    Parameters
    ---------- 
    datain : array_like 
        A 2D data array (IRIS slit)
    psfin  : array_like 
        The point spread function to be applied in the y-direction
    div : bool 
        Set to True to divide in Fourier space (Default is False, so multiply in Fourier space)
    rev_psf : bool 
        Set to reverse the 1D input PSF
              
    Returns
    -------            
    dataout : array_like 
        The input data convolved with the PSF
    
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
        pin = pin/(np.sum(pin))
        
    #Pad the PSF if it is too short
    if ydiff > 0:
        rs = ydiff//2
        pin=np.pad(psfin, [rs, rs+ydiff%2], 'edge')
        #Use edge values or you get a divide by zero error
        #renormalize PSF (dx=1)
        pin = pin/(np.sum(pin))

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
