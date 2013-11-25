
pro seg_rowval, infile, segment, row, value, FULL=full, $
  MIN_COL=min_col, MAX_COL=max_col, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

;+
; NAME:
;   SEG_ROWVAL
;
; PURPOSE:
;   Set the value of a single row in a FITS file extension. Optionally
;   include only a certain range of columns in the segment. 
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    SEG_ROWVAL, infile, segment, col, value, [ VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    COL: column number to set
;    VALUE: value to set the pixels to
;
; OPTIONAL INPUTS:
;    NONE
;
; INPUT KEYWORD PARAMETERS:
;    FULL : if set, do the full row, including pre and over scan areas
;    MIN_COL : minimum column number at which to set the value
;    MAX_COL : maximum column number at which to set the value
;    PRE_COLS:  number of prescan columns in image segment, default=10
;    OVER_COLS: number of overscan columns in image segment, default=20
;    OVER_ROWS: number of overscan rows in image segment, default=20
;    VERBOSE: if set, print some diagnostics
;    
; OUTPUT KEYWORD PARAMETERS:
;    NONE
; 
; HISTORY:
; Written by: P. Doherty  April, 2013
;
;-
;-----------------------------------------------------------------------------
compile_opt IDL2

IF N_Elements(pre_cols)  EQ 0 THEN pre_cols = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

; check for a readable file and legitimate segment number
errmsg='No error.'
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

; read in the image data from the input file
fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; check for legitimate row number request
if ((row GT s[2]-1) || (row LT 0)) then begin
  print,'SEG_ROWVAL: Invalid row number.'
  return
endif

; by default, only allow work in the imaging area
IF N_Elements(min_col) EQ 0 THEN min_col = pre_cols
IF N_Elements(max_col) EQ 0 then max_col = s[1]-(over_cols+1)

; do full row if requested, including serial prescan and overscan
; ignoring MIN_COL and MAX_COL values
if keyword_set(FULL) then begin
  min_col = 0
  max_col = s[1]-1
endif

; check to be sure the col numbers are acceptable
if ((min_col LT 0) || (min_col GT s[1]-1)) then begin
  print,'SEG_ROWVAL: Invalid column number.'
  return
endif

if ((max_col LT 0) || (max_col GT s[1]-1)) then begin
  print,'SEG_ROWVAL: Invalid column number.'
  return
endif

; make the change and save the file
im[min_col:max_col,row] = value

if file_test(infile, /WRITE) then MODFITS, infile, im, hdr, EXTEN_NO=segment

end

