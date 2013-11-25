
pro seg_colval, infile, segment, col, value, FULL=full, $
  MIN_ROW=min_row, MAX_ROW=max_row, $
  PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
  VERBOSE=verbose

;+
; NAME:
;   SEG_COLVAL
;
; PURPOSE:
;   Set the value of a column in a FITS file extension. 
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    SEG_COLVAL, infile, segment, col, value, [ VERBOSE=verbose ]
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
;    FULL : if set, do the whole column, including overscan pixels
;    MIN_ROW : the lowest row in the column to set
;    MAX_ROW : the highest row in the column to set
;    PRE_COLS : number of prescan columns pre row 
;    OVER_COLS : number of overscan columns present
;    OVER_ROWS : number of overscan rows present 
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

IF N_Elements(min_row)   EQ 0 THEN min_row = 0
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

if keyword_set(VERBOSE) then verbose = 1 else verbose = 0

; check that the file is readable, segment number OK, and such
errmsg='No error.'
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

; read in the input file
fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; check for a legitimate column number request
if ((col GT s[2]-1) || (col LT 0)) then begin
  print,'SEG_COLVAL: Invalid column number.'
  return
endif

; if max_row not specified, set it to the top of imaging area or segment
IF N_Elements(max_row) EQ 0 then begin
    ; by default, only allow work in the imaging area, not parallel overscan
    max_row = s[2]-(over_rows+1)
    ; in prescan or overscan do full column
    if ((col LT pre_cols) || (col GE s[1]-(over_cols))) then max_row = s[2]-1 
endif

; Do full column if requested, including parallel overscan, 
; ignoring MIN_ROW and MAX_ROW values
if keyword_set(FULL) then begin
   min_row = 0
   max_row = s[2]-1 
endif

; check to be sure the row numbers are acceptable
if ((min_row LT 0) || (min_row GT s[2]-1)) then begin
  print,'SEG_COLVAL: Invalid row number.'
  return
endif
if ((max_row LT 0) || (max_row GT s[2]-1)) then begin
  print,'SEG_COLVAL: Invalid row number.'
  return
endif

; make the change and save the file
im[col, min_row:max_row] = value

if file_test(infile, /WRITE) then MODFITS, infile, im, hdr, EXTEN_NO=segment

end

