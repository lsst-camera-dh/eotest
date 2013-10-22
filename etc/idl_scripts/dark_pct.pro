

; dark_pct return the max dark current for a given percentile of pixels
function dark_pct, file1, file2, file3, percentile, GAIN=gain, EXT=ext

compile_opt defint32, strictarr

if ~keyword_set(EXT)  then ext=1
if ~keyword_set(GAIN) then gain=1.0
if ~keyword_set(COUNT) then count = 1000

percentile = float(percentile)

if ((percentile GT 100.0) || (percentile LT 0)) then begin
    print, "Error. Percentile must be > 0 and < 100."
    return, -1
endif

im = float(fits_median(file1, file2, file3))
if(n_elements(im) EQ 1) then begin
   print,"Error. Image median process failed"
   return, -1
endif

; get exposure time from one of the files
header = HEADFITS(file1 , /SILENT)
exptime = SXPAR( header, 'EXPTIME', /SILENT)

bmean = mean(float(im[525:530, 100:2000]))    ; get bias from overscan
im = im[10:521, 0: 2002]               ;  the part of the image we want to test
im = ((im - bmean) * gain) / exptime   ; electrons per pixel per second
pc  = im[ (SORT(im))[percentile * N_ELEMENTS(im) / 100] ]

return, pc
end
