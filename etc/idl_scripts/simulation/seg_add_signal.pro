
; poisson noise is disabled until I figure out how to do it correctly

pro seg_add_signal, infile, segment, signal, GAIN=gain, POISSON=poisson, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

;+
; NAME:
;   SEG_ADD SIGNAL
;
; PURPOSE:
;   Add signal to an image segment. 
;
; EXPLANATION:
;   Add signal to an image segment. Signal is defined in electrons and
;   an optional gain parameter may be used for scaling.
;   
; CALLING SEQUENCE:
;    SEG_ADD_SIGNAL, infile, segment, signal [ , GAIN=gain, POISSON=poisson,
;       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
;       VERBOSE=verbose
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    SIGNAL: signal, in electrons, to add
;
; OPTIONAL INPUTS:
;    GAIN      : system gain in e-/DN
;    POISSON   : if set, add poisson noise to the signal
;    PRE_COLS  : number of prescan columns in image segment, default=10
;    OVER_COLS : number of overscan columns in image segment, default=20
;    OVER_ROWS : number of overscan rows in image segment, default=20
;    VERBOSE   : if set, print some diagnostics
;
; INPUT KEYWORD PARAMETERS:
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

IF N_Elements(gain)      EQ 0 THEN gain = 1.0      ; system gain in e-/DN
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20
if keyword_set(poisson) then poisson = 1 else poisson = 0
if keyword_set(verbose) then verbose = 1 else verbose = 0
  
; disable poisson for now until I figure it out
poisson=0

errmsg=' '
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

signal = float(signal)

fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; create an image array to hold the calculated signal pixels
; (just the image aray pixels)
sigim = float(im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] )

; zero the signal image
sigim[*] = 0.00

; set the signal image pixels
sigim[*] = float(signal)

; disable poisson for now until I figure it out
poisson=0
if keyword_set(poisson) then begin
   ; add noise using the 'noise' function with the /poisson keyword set
   ; sigim = noise(sigim,/poisson)  THIS SOMETIMES CRASHES!
   
   sigim = sigim + sqrt(mean(sigim)) * PoiDev(sigim)
endif
   

sigim = sigim/gain

; add signal image into image array of input image
im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] = $
    im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))]  + sigim

im = uint(im)

; save the file
MODFITS, infile, im, hdr, EXTEN_NO=segment

end

