
pro seg_noise, infile, segment, noise, GAIN=gain, $
  FULL=full, IMAGE=image, SER_PRE=ser_pre, SER_OVER=ser_over, PAR_OVER=par_over, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

;+
; NAME:
;   SEG_NOISE
;
; PURPOSE:
;   Add noise to the pixels in a FITS file extension. Optionally works on 
;   image area, parallel overscan, serial prescan, and serial overscan
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    SEG_NOISE, infile, segment, noise [ , GAIN=gain, $
;      FULL=full, IMAGE=image, SER_PRE=ser_pre, SER_OVER=ser_over, PAR_OVER=par_over, $
;      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
;      VERBOSE=verbose
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    NOISE: rms value of the noise distribution
;
; OPTIONAL INPUTS:
;
; INPUT KEYWORD PARAMETERS:
;    FULL      : if set, add noise to all the pixels
;    IMAGE     : if set, add noise to pixels in imaging array
;    SER_OVER  : if set, add noise to pixels in serial overscan columns
;    SER_PRE   : if set, add noise to pixels in serial prescan columns
;    PAR_OVER  : if set, add noise to pixels in parallel overscan rows
;    PRE_COLS  : number of prescan columns in image segment, default=10
;    OVER_COLS : number of overscan columns in image segment, default=20
;    OVER_ROWS : number of overscan rows in image segment, default=20
;    VERBOSE   : if set, print some diagnostics
;
;
; OUTPUT KEYWORD PARAMETERS:
;
; EXAMPLES:
;
; WARNING:
;
; PROCEDURES USED:
;
; HISTORY:
; Written by: P. Doherty  April, 2013
;
;-
;-----------------------------------------------------------------------------
compile_opt IDL2

IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20
IF N_Elements(gain) EQ 0 THEN gain = 1.0 ; in e-/DN
 
errmsg=' '
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

gain  = double(gain)
noise = double(noise)

fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; if PAR_OVER then add noise to the parallel overscan region 
; Note: parallel overscan does not include serial prescan or overscan areas
if keyword_set(PAR_OVER) then begin
   im2 = double(im[pre_cols:(s[1]-(over_cols+1)), (s[2]-over_rows):(s[2]-1)] )
   s2 = size(im2)
   im3 = randomn(seed,s[1], s[2], /double) * (double(noise)/double(gain))
   im[pre_cols:(s[1]-(over_cols+1)), (s[2]-over_rows):(s[2]-1)]  = uint(im2 + im3)
endif

; if SER_OVER then add noise to the serial overscan
if keyword_set(SER_OVER) then begin
   im2 = double(im[(s[1]-over_cols):(s[1]-1), 0:(s[2]-1)])
   s2 = size(im2)
   im3 = randomn(seed,s[1], s[2], /double) * (double(noise)/double(gain))
   im[(s[1]-over_cols):(s[1]-1), 0:(s[2]-1)] = im2 + im3
endif
  
; if SER_PRE then add noise to the serial prescan
if keyword_set(SER_PRE) then begin
   im2 = double(im[0:(pre_cols-1), 0:s[2]-1])
   s2 = size(im2)
   im3 = randomn(seed,s[1], s[2], /double) * (double(noise)/double(gain))
   im[0:(pre_cols-1), 0:s[2]-1] = im2 + im3
endif

; if IMAGE then add noise to the imaging array
if keyword_set(IMAGE) then begin
   im2 = double(im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))])
   s2 = size(im2)
   im3 = randomn(seed,s[1], s[2], /double) * (double(noise)/double(gain))
   im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] = uint(im2 + im3)
endif

; if FULL is set, do the whole segment
if keyword_set(FULL) then begin
   im2 = double(im)
   s2 = size(im2)
   im3 = randomn(seed,s[1], s[2], /double) * (double(noise)/double(gain))
   im = uint(im2 + im3)
endif

; save the file
MODFITS, infile, im, hdr, EXTEN_NO=segment

end




