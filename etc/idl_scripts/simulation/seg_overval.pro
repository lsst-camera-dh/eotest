
pro seg_overval, infile, segment, value, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

;+
; NAME:
;   SEG_OVERVAL
;
; PURPOSE:
;   Set the value of all prescan and overscan area pixels in a FITS 
;   file extension. 
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    SEG_OVERVAL, infile, segment, value, [ PRE_COLS=pre_cols, 
;        OVER_COLS=over_cols, OVER_ROWS=over_rows, VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    VALUE: value to set the pixels to
;
; OPTIONAL INPUTS:
;    NONE
;
; INPUT KEYWORD PARAMETERS:
;    PRE_COLS:  number of prescan columns in image segment, default=10
;    OVER_COLS: number of overscan columns in image segment, default=20
;    OVER_ROWS: number of overscan rows in image segment, default=20
;    VERBOSE: if set, print some diagnostics
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
  
nParam = N_Params()
if nParam NE 3 then begin
  print,'SEG_OVERSCAN: Incorrect number of parameters. Exiting.'
  return
endif

errmsg=' '
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

seg_imval,infile,segment,value,/SER_PRE,/SER_OVER,/PAR_OVER

end

