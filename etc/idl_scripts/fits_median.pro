; fits_median : return an image which is the median of three images

function fits_median, file1, file2, file3, FIX=fix

if ~keyword_set(FIX) then fix = 0

header1 = HEADFITS(file1 , /SILENT)  ; get master header for file1
header2 = HEADFITS(file2 , /SILENT)  ; get master header for file2
header3 = HEADFITS(file3 , /SILENT)  ; get master header for file3

time1 = SXPAR( header1, 'EXPTIME', /SILENT)
time2 = SXPAR( header2, 'EXPTIME', /SILENT)
time3 = SXPAR( header3, 'EXPTIME', /SILENT)

if (time1 NE time2) || (time2 NE time3) ||  (time1 NE time3) then begin
   print,"Error: unequal exposure times. "
   return, -1
endif 

fits_read, file1, im1, hdr1 ;, EXTEN_NO=ext
fits_read, file2, im2, hdr2 ;,; EXTEN_NO=ext
fits_read, file3, im3, hdr3  ;,; EXTEN_NO=ext

header1 = HEADFITS(file1, /SILENT)
header2 = HEADFITS(file2, /SILENT)
header3 = HEADFITS(file3, /SILENT)

if (fix EQ 1) then begin
   m1 = median(im1)
   m2 = median(im2)   
   m3 = median(im3)
   med = (m1 + m2 + m3 ) / 3.0
   err1 = m1 - med
   err2 = m2 - med
   err3 = m3 - med
   im1 = im1 - err1
   im2 = im2 - err2
   im3 = im3 - err3
endif

; Generate image that is median of the three input images
imgarray = [ [[im1]], [[im2]], [[im3]] ]               ; create 3d array with the images
totalImage = Total(imgarray, 3)                        ; total the images in the dim 3.
minValue = Min(imgarray, Dimension=3)                  ; 2d array of min values
maxValue = Max(imgarray, Dimension=3)                  ; 2 d array of max values
medianImage = totalImage - minValue - maxValue         ; med- total - min- max

return, medianImage

end