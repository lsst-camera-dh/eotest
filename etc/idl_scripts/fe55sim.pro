

pro fe55sim, outfile, XSIZE=xsize, YSIZE=ysize, EVENTS=events, BETAS=betas, $
             GAIN=gain, OFFSET=offset, NOISE=noise

if ~keyword_set(XSIZE)  then xsize  = 542     ; # of columns in output image
if ~keyword_set(YSIZE)  then ysize  = 2022    ; y size of output image
if ~keyword_set(EVENTS) then events = 10000   ; number of x-ray events
if ~keyword_set(BETAS)  then betas  = 0       ; fraction of beta events to include (~12% in real world)
if ~keyword_set(GAIN)   then gain   = 1.00    ; system gain e-/DN
if ~keyword_set(OFFSET) then offset = 0       ; offset ot add to the image after making events
if ~keyword_set(NOISE)  then noise  = 0       ; peak to peak noise

if (xsize GT 10000) or  (xsize LT 1) then begin
    print, 'xsize must be > 1  and  < 10000'
    exit
endif

if (ysize GT 10000) or  (ysize LT 1) then begin
    print, 'ysize must be > 1  and  < 10000'
    exit
endif

if (betas GT 1) or  (betas LT 0) then begin
    print, 'beta fraction must be > 0  and  < 1'
    exit
endif

alphas = events * (1 - betas)
betas  = events * (betas)
print,alphas,betas

im=uintarr(xsize, ysize)
im[*] = offset

s=size(im)
print, s

if (noise GT 0) then begin
   ; generate a random noise image
   s = size(im)
   seed = systime(/seconds)
   noise_im = float(randomn(seed, s[1], s[2]));
   ; make the noise image have the correct peak to peak range
   range = max(noise_im) - min(noise_im)
   noise_im = noise_im * (noise / range)
   ; make the noise image have a mean of zero
   offset = mean(noise_im)
   noise_im = noise_im - offset
   ; generate the combined data+noise image
   im = float(im) + noise_im
endif

wait, 1 & seed = systime(/seconds) MOD 1000 ; reasonably random seed

if (alphas GE 1) then begin
   print, 'Generating ', alphas, ' alpha events'
   for i = 1, alphas do begin
      x = uint(randomu(seed) * (xsize - 20))    ; generate random event location
      y = uint(randomu(seed) * (ysize - 20))
      im[x, y] = im[x, y] + (1620 / gain)
   endfor
endif

if (betas GE 1) then begin
   print, 'Generating ', betas, ' beta events'
   for i = 1, betas do begin
      x = uint(randomu(seed) * (xsize - 20))    ; generate random event location
      y = uint(randomu(seed) * (ysize - 20))
      im[x, y] = im[x, y] + (1778 / gain)
   endfor
endif

im[where(im LT 0)] = 0
im = uint(im)
fits_write, outfile, im

end




            