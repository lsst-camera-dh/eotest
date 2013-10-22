
; ctesim : degrade an image by simulating CTE at a specified level
pro ctesim, infile, outfile, PCTE=pcte, SCTE=scte

; should verify that the parameters are all correct here
if ~keyword_set(PCTE) then pcte=1.00
if ~keyword_set(SCTE) then scte=1.00

pcte = float(pcte)
pcti = 1.00 - pcte

scte = float(scte)
scti = 1.00 - scte

fits_read, infile, im1, hdr
im1=float(im1)
s=size(im1)
cols = s[1]
rows = s[2]
im2 = fltarr(cols,rows)
im2[*,*] = 0.00
outim = fltarr(cols,rows)
outim[*,*] = 0.00

; subtract off bias median before processing
med = median(im1[530:540,0:2000])
im1 = im1 - med

outim=im1   ; may or may not get processed, so set output file to input file

if (PCTE NE 1.00) then begin
   for i = 0, rows-2 do begin
      if (( i MOD 100) EQ 0 ) then print, "processing row i = ", i
      outim[*, i] = im1[*,0] * pcte            ; copy bottom row to output
      for j = 0, rows-2 do begin               ; then calculate a new shifted frame
        im2[*, j] = (im1[*, j] * pcti) + (im1[*,j+1] * pcte)
      endfor
      im2[*, rows-1] = (im1[*, rows-1] * pcti) ; do the last row a little differently
      im1 = im2                                ; done making new image, copy it back to im1
      im1[where(im1 LT 0.0001)] = 0
   endfor
endif

if (SCTE NE 1.00) then begin
   im1 = outim
   for i = 0, cols-2 do begin
      if (( i MOD 100) EQ 0 ) then print, "processing column i = ", i
      outim[i, *] = im1[0,*] * scte  ; copy column to output
      ; then calculate a new shifted frame
      for j = 0, cols-2 do begin
        im2[j, *] = (im1[j, *] * scti) + (im1[j+1,*] * scte)
      endfor
      
      im2[cols-1, *] = (im1[cols-1, *] * scti) ; do the last row a little differently
      im1 = im2                                ; done making new image, copy it back to im1
      im1[where(im1 LT 0.0001)] = 0            ; clip off annoying small fractions
   endfor
endif

outim = outim + med          ; add median back in
outim = uint(outim)          ; convert to integers (electrons are quantized, no?)
fits_write, outfile, outim   ; write the output file
print,'Done.'
end