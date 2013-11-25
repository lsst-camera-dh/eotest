
pro seg_imval, infile, segment, value, $
  FULL=full, IMAGE=image, SER_PRE=ser_pre, SER_OVER=ser_over, PAR_OVER=par_over, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

;+
; NAME:
;   SEG_IMVAL
;
; PURPOSE:
;   Set the value of pixels in a FITS file extension. Optionally works on 
;   image area, parallel overscan, serial prescan, and serial overscan
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    SEG_IMVAL, infile, segment, value, [ SER_OVER=ser_over, PAR_OVER=par_over, 
;        PRE_COLS=pre_cols, OVER_COLS=over_cols,OVER_ROWS=over_rows, 
;        VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    VALUE: value to set the pixels to
;
; OPTIONAL INPUTS:
;    IMAGE     : if set, set values of pixels in imaging array
;    SER_OVER  : if set, set values of pixels in serial overscan columns
;    SER_PRE   : if set, set values of pixels in serial prescan columns
;    PAR_OVER  : if set, set values of pixels in parallel overscan rows
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

IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20
  
errmsg=' '
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; if PAR_OVER then set the value of the parallel overscan region 
; Note: parallel overscan does not include serial prescan or overscan areas
if keyword_set(PAR_OVER) then $
   im[pre_cols:(s[1]-(over_cols+1)), (s[2]-over_rows):(s[2]-1)] = value

; if SER_OVER then set the value of the serial overscan
if keyword_set(SER_OVER) then $
   im[(s[1]-over_cols):(s[1]-1), 0:(s[2]-1)] = value
  
; if SER_PRE then set the value of the serial prescan
if keyword_set(SER_PRE) then $
   im[0:(pre_cols-1), 0:s[2]-1] = value

; if IMAGE then set the value of the imaging array
if keyword_set(IMAGE) then $
   im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] = value

; if FULL is set, do the whole segment
if keyword_set(FULL) then $
   im[0:s[1]-1, 0:s[2]-1] = value
 
; save the file
MODFITS, infile, im, hdr, EXTEN_NO=segment

end

