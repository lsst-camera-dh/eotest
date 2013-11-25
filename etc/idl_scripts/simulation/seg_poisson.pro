
pro seg_poisson, infile, segment, GAIN=gain, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

;+
; NAME:
;   SEG_POISSON
;
; PURPOSE:
;   Add poisson noise to the imaging array pixels in a FITS file extension. 
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    SEG_POISSON, infile, segment,  [ PRE_COLS=pre_cols, OVER_COLS=over_cols,
;        OVER_ROWS=over_rows, VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    NOISE: rms value of the noise distribution
;
; OPTIONAL INPUTS:
;
; INPUT KEYWORD PARAMETERS:
;    IMAGE     : if set, add noise to  pixels in imaging array
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

fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; extract the image aray pixels, subtract the overscan bias offset
; add noise to the result using the 'noise' function with the /poisson
; keyword set, add the offset back in, put the pixels back into the
; imaging array

;extract the image aray pixels
im2 = im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] 

;subtract the overscan bias offset
offset = seg_mean(infile, segment, /ser_over)
im3 = im2 - offset ; adjust for offset
im3 = im3 / gain   ; adjust for gain

; add noise using the 'noise' function with the /poisson keyword set
im3 = noise(im3,/poisson)

; take out the DC component, we want just the noise
im3 = im3 - mean(im3)

; put the noise into the pixels in the imaging array
im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] = im2 + im3

; save the file
MODFITS, infile, im, hdr, EXTEN_NO=segment

end




