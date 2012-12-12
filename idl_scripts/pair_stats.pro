
; pair_stats:
; calculate the gain and noise of a CCD camera system by examining
; two flat field and two bias frames. The calulation is the standard
; mean/variance thing.

; result is an array of floats:
; bias mean, bias std dev, flat mean, flat variance, gain, noise

pro pair_stats, flat, ext, result
    ;print, flat, ext
    infile1 = string(flat,'_flat1.fits', format='(A,A)')
    infile2 = string(flat,'_flat2.fits', format='(A,A)')
    ;bfile1 = string(bias,'_bias1.fits', format='(A,A)')
    ;bfile2 = string(bias,'_bias2.fits', format='(A,A)')
      
    ; read in the fits files
    fits_read, infile1, flat1, fhdr1, EXTEN_NO=ext
    fits_read, infile2, flat2, fhdr2, EXTEN_NO=ext
    ;fits_read, bfile1, bias1, bhdr1, EXTEN_NO=ext
    ;fits_read, bfile2, bias2, bhdr2, EXTEN_NO=ext
    ; use image area as bias regions for computations
    ;b1 = float(bias1[200:300, 900:1000])
    ;b2 = float(bias2[200:300, 900:1000])
    ; use overscan as bias regions for computations
    b1 = float(flat1[525:530, 100:2000])
    b2 = float(flat2[525:530, 100:2000])
    ; calculate the bias mean
    bmean = (mean(b1) + mean(b2)) / 2  
    ; extract the desired flat field regions for computations   
    f1 = float(flat1[200:300, 900:1000] - bmean)
    f2 = float(flat2[200:300, 900:1000] - bmean)
    ; calculate the flat ratio and correct for it
    fratio = float(mean(float(f1)/float(f2)))   
    f2 = float(f2) * fratio       
    ; calculate the mean value of the flat field images
    fmean = (mean(f1) + mean(f2)) / 2  
    ; create flat field and bias difference images
    fdiff = f1 - f2               
    bdiff = b1 - b2
    ; calculate the variances (std dev of diff image squared over 2)
    fvar = (stddev(fdiff)^2)/2  
    bvar = (stddev(bdiff)^2)/2
    ; calculate gain in e-/DN
    gain =  fmean / fvar    
    ; calculate noise in e- rms    
    bias_rms = stddev(b1)
    noise =  float(gain * (stddev(b1)))  
    
    ; the results to be returned
    result = [bmean, bias_rms, fmean, fvar, gain, noise]
    ;print,flat,result
end
