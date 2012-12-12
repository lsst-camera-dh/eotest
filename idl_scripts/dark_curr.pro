

; dark_curr return the dark current based on three images.
function dark_curr, file1, file2, file3, EXT=ext, GAIN=gain, COUNT=count

if ~keyword_set(EXT)  then ext=1
if ~keyword_set(GAIN) then gain=1.0
if ~keyword_set(COUNT) then count = 1000

dcurrent = fltarr(count)

im = fits_median(file1, file2, file3)
if(n_elements(im) EQ 1) then begin
   print,"Error. Image median process failed"
   return, -1
endif
   
bmean = mean(float(im[525:530, 100:2000]))    ; get bias from overscan

; perform measurement in many random locations, then take median result
; do this to avoid measuring where bright defects exist
wait, 1 & seed = systime(/seconds) MOD 1000 ; reasonably random seed
for i = 1, count do begin
   x = uint(randomu(seed) * 400)    ; generate random box location
   y = uint(randomu(seed) * 1900)
   signal = mean(float(im[x+10:x+100, y:y+100] - bmean))  ; get array signal level
   dcurrent[i-1] =  signal 
endfor
dark_current = median(dcurrent) * gain

return, dark_current
end
