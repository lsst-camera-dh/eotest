function xray_gain, file, EXT=ext, BINSIZE=binsize, REGION=region
compile_opt idl2
 On_error,2

 if ~keyword_set(EXT ) then ext = 1
 if ~keyword_set(BINSIZE) then binsize = 10
 if ~keyword_set(REGION) then region = 9

 ; check if the file exists and is readable , if not exit
 if (~file_test(file, /READ)) then begin
    print,'Error: Specified file does not exist or is unreadable'
    return, -1
 endif

 fits_read,file,image,hdr,exten_no=ext
 image=float(image)
 noise = stddev(image[535:540,200:900]) ; note hard coded region!
 ; pixels are considered part of an x-ray event if they have a value greater
 ; than the image median plus some margin. Cycle through the calculation several
 ; times using 'margin' values that are multiples of the image noise. At the
 ; end, take the median of the computed gain values and erport that as the result
 gains= fltarr(10)
 for j = 0, 9 do begin
    ; compute margin value and threshold for identifying x-ray events
   margin = noise * (j+3)
   threshold = median(image) + margin
   ; create a binary image by thresholding the x-ray image then label each
   ; contiguous region in the image and make a histogram that indicates how
   ; many pixels are in each region by counting the pixels with that region
   ; index. Save the reverse indices to use later in extracting pixel values
   ; from the original image
   event_im = image GT threshold
   label_im = label_region(event_im)
   region_sizes=histogram(label_im, min=0, reverse_indices=R)
   region_count = n_elements(region_sizes)
   ; Subtract the median pixel value from the image. Then create an array of
   ; summed pixel values for each region using the reverse indices to reach
   ; back into the image. The computed 'signal' is an array of values
   ; each represeting the total charge in a single x-ray event.
   image = image - median(image)
   signal = fltarr(region_count)
   for i = 0, region_count-1 do begin
      p = R[R[i]:R[i+1]-1]
      signal[i] = total(image[p])
   endfor
   ; Extract the regions that are smaller than the minimum size specified via
   ; the region keyword. If no regions are found, then set this measurement to
   ; 0 so it will not affect the result
   event_index = where(region_sizes LE region)
   if (event_index[0] EQ -1) then begin
      print,'Error: No x-ray events found'
      gains[j] = 0
      break
   endif
   ; Calculate a histogram of the signal values, and find peak value location
   ; Compute the DN value associated with peak location accounting for the
   ; histogram bin size and the margin
   spectrum = histogram(signal[event_index], min=0, binsize=binsize)
   peak = mean(where(spectrum EQ max(spectrum)))
   signal= float(peak) * float(binsize) + float(margin)
   ; Compute system gain by dividing known number of electrons per event by
   ; measured DN per event
   gains[j] = 1620.0/signal
endfor
; compute median result value and return it
gain = median(gains)
return, gain
end
