pro seg_expose, infile, segment, flux, exptime, $
   GAIN=gain, POISSON=poisson, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
; NAME:
;   SEG_EXPOSE
;
; PURPOSE:
;   Add simulated flux to an image segment
;
; EXPLANATION:
;   Simulate an exposure of flux for a certain time, by adding
;   signal to an image segment. Flux is in electrons per second. Time
;   is in seconds.
;   
; CALLING SEQUENCE:
;    SEG_EXPOSE, infile, segment, flux, exptime, [ , GAIN=gain, 
;       POISSON=poisson, PRE_COLS=pre_cols, OVER_COLS=over_cols, 
;       OVER_ROWS=over_rows, VERBOSE=verbose
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    FLUX: electrons per second
;    EXPTIME: exposure duration in seconds
;
; OPTIONAL INPUTS:
;    GAIN      : system gain in e-/DN
;    POISSON   : if set, ad poisson noise to the signal
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
 
if keyword_set(POISSON) then poisson = 1 else poisson = 0
if keyword_set(VERBOSE) then verbose = 1 else verbose = 0

if N_Params() NE 4 then begin
  print,'SEG_EXPOSE: Incorrect number of parameters. Exiting.'
  return
endif

errmsg=' '
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

signal = float(flux) * float(exptime)

; update the header with useful information
fits_read, infile, im, hdr, exten_no=segment
SXADDPAR, seghdr, 'GAIN    ', gain, 'segment gain',  FORMAT='(F8.3)'
MODFITS, infile, 0, hdr, EXTEN_NO=segment

; add the signal to the image extension
seg_add_signal, infile, segment, signal, GAIN=gain, POISSON=poisson, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

end

