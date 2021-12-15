import numpy as np
import scipy.fft as fft


def IRIS_SG_deconvolve(data_in, psf, 
                       iterations=10, 
                       fft_div=False):


    '''

    Graham S. Kerr
    July 2020
 
    NAME:   IRIS_SG_Deconvolve.py

    PURPOSE: Deconvolves IRIS SG data using the PSFs from Courrier et al 2018.
    
    INPUTS: data_in -- A 2D IRIS SG array [ypos, wavelength]
            psf     -- The appropriate PSF 
                       These are not currently in iris_lmsalpy, so I have just 
                       saved the IDL versions and restore them in the notebook 
                       before I call this function. 
            iterations -- The number of Richardson Lucy iterations to run through
                          Default = 10
            fft_div -- Set to use skip iterations and instead deconvolve by division
                       in Fourier Space   

    NOTES: Based on iris_sg_deconvolve.pro by Hans Courrier, but not all the functionality
           is included here yet

           There are probably more clever ways to code this in -- i'm fairly new to python. 

            To Do: Add error statements

    '''


    # Remove negative values 
    dcvim = data_in.copy()
    dcvim[np.where(dcvim<0)] = 0
    data_in_zr = dcvim

    if fft_div == True:
        dcvim = FFT_conv_1D(data_in,psf,div=True)
    else:
        for ind in range(1,iterations+1):
            #print('iteration = %3d' %(ind))
            step1 = data_in_zr/(FFT_conv_1D(dcvim,psf,rev_psf=False,div=False))
            #print(np.nanmax(step1[265,:]))
            step2 = FFT_conv_1D(step1,psf, rev_psf=True)
            dcvim = dcvim * step2

    return dcvim


def FFT_conv_1D(datain, psfin, div = False, rev_psf=False):
    

    '''

    Graham Kerr
    July 2020
   
    NAME: FFT_conv_1D
    
    PURPOSE: Function to do FFT convolution in the y-direction 
              of the input data (first dimensioin). This way we 
              can pass the 1D PSF.
          
    INPUTS: datain -- a 2D data array [nominally, slit pos vs wavelength]
            psfin  -- the PSF to be applied in the y-direction
            imsize -- the dimensions of the input data array 
            psflen -- the length of the psf. 

    KEYWORDS: div -- Set to True to divide in Fourier space. 
                     Default is False, so multiply in Fourier space.
              rev_psf -- Set to reverse the 1D input PSF
              
            
    OUTPUTS: dataout -- the input data convolved with the PSF
    
    NOTES: Pretty much copied exactly from Hans Courrier's IDL version 
            in the SSW IRIS software tree, as part of iris_sg_deconvole.pro
           
            Can probably be written in a much better way more suitable for 
            python.

    '''
  
    # length of input psf
    psflen = len(psfin)
    
    # dimensions of input data
    imsize = datain.shape
    
    # Get difference of image size and psf length
    ydiff = imsize[0]-psflen  
    
    # Cut the PSF if it is too long
    if ydiff <= 0:

        rs = int(np.abs(ydiff)/2)

        if np.abs(ydiff) % 2 == 1:
            pin = psfin[rs+1:psflen-rs] # odd ydiff
        else:
            pin =  psfin[rs:psflen-rs] # even ydiff

       # renormalize PSF
        pin = pin/np.sum(pin) 
        
    # Pad the PSF if it is too short
    if ydiff > 0:
        
        rs = int(ydiff/2)
        
        padl = np.zeros(rs,dtype=float)
        
        if ydiff % 2 == 1:
            padr = np.zeros(rs+1,dtype=float)
        else:
            padr = np.zeros(rs,dtype=float)
            
        pin = np.concatenate((padl,psfin,padr))
       
    # Replicate the PSF over wavelength array also
    pin_full = np.tile(pin,[imsize[1],1])
    pin_full = np.transpose(pin_full)
    
    # Shift PSF center to zero to center output
    pin_full = np.roll(pin_full, (int(-imsize[0]/2),0),axis=0)
        
    # Reverse the PSF if needed
    if rev_psf == True:
        pin_full = np.flip(pin_full,axis=0)        
        if psflen % 2 == 0:
            pin_full = np.roll(pin_full,(1,0),axis=0)
            
    
    # Perform the FFT
    fpsf = fft.fft(pin_full,axis=0)
    datain=datain.astype(np.float64)
    fdatain = fft.fft(datain,axis=0)
    
    
    # Multiply(divide) the PSF and data, and convert back to k space
    if div == False: 
        dataout = fft.ifft((fdatain*fpsf),axis=0).real
    else:
        dataout = fft.ifft((fdatain/fpsf),axis=0).real
    
    
    return dataout
