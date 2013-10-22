function flat_gain, file1, file2, EXT=ext, COUNT=count
;+
; NAME:
;       FLAT_GAIN
; PURPOSE:
;         
; EXPLANATION:
; calculate the gain and noise of a CCD camera system by examining
; two flat field imagees. The calulation is the standard
; mean/variance thing.
   
;
; CALLING SEQUENCE:
;       
;
; INPUT PARAMETERS:
;       
;
; OPTIONAL OUTPUTS:
;       None
;
; OPTIONAL KEYWORD INPUTS:
;       
; OPTIONAL KEYWORD OUTPUT:
;       
;
; PROCEDURE:
;       
;
; PROCEDURE CALLS:
;       
; MODIFICATION HISTORY:
;       
;- 
;-------------------------------------------------------------------------------

if ~keyword_set(EXT) then ext = 1         ; image extension to process
if ~keyword_set(COUNT) then count = 1000  ; number of measurements to make

; check if the files exists and is readable , if not exit
if (~file_test(file1, /READ)) then begin
     print,'Error: Specified file does not exist or is unreadable'
     return, -1
endif
if (~file_test(file2, /READ)) then begin
    print,'Error: Specified file does not exist or is unreadable'
    return, -1
endif

; read in the fits files
fits_read, file1, flat1, fhdr1, EXTEN_NO=ext
fits_read, file2, flat2, fhdr2, EXTEN_NO=ext

; extract the desired bias regions for computations, calculate mean value
b1 = float(flat1[535:540, 100:2000])
b2 = float(flat2[535:540, 100:2000])
bmean = (mean(b1) + mean(b2)) / 2 

gains=fltarr(count)
seed = systime(/seconds) MOD 1000 ; reasonably random seed

for i = 1, count do begin
   ; generate random box location
   x = uint(randomu(seed) * 400)
   y = uint(randomu(seed) * 1900)
   
   ; extract the desired flat field regions for computations   
   f1 = float(flat1[x:x+100, y:y+100] - bmean)
   f2 = float(flat2[x:x+100, y:y+100] - bmean)
   ; calculate the flat ratio and correct for it
   fratio = float(mean(float(f1)/float(f2)))   
   f2 = float(f2) * fratio       
   ; calculate the mean value of the flat field images
   fmean = (mean(f1) + mean(f2)) / 2  
   ; create flat field and bias difference images
   fdiff = f1 - f2               
   ; calulate the variance of the flat difference image
   fvar = (stddev(fdiff)^2)/2  
   ; calculate gain in e-/DN
   gains[i-1] =  fmean / fvar     

endfor
gain = median(gains)
return, gain

end



