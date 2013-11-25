
pro seg_pixval, infile, segment, col, row, value, VERBOSE=verbose

;+
; NAME:
;   SEG_PIXVAL
;
; PURPOSE:
;   Set the value of pixels in a FITS file extension. 
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    SEG_PIXVAL, infile, segment, col, row, value, [ VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    COL: column number of the pixel to set, can be an array
;    ROW: row number of the pixel to set, can be an array
;    VALUE: value to set the pixels to
;
; OPTIONAL INPUTS:
;    NONE
;
; INPUT KEYWORD PARAMETERS:
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

errmsg='No error.'
if (seg_valid(infile, segment, /write,  message=errmsg) NE 0) then begin
  print,errmsg
  return
end

ncol = n_elements(col)
nrow = n_elements(row)

if (n_elements(col) NE n_elements(row)) then begin
   print,'SEG_PIXVAL: Mismatched row and column arrays. Exiting.'
   return
endif

fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; check for valid row and column numbers
if ((max(col) GT s[1]-1) || (min(col) LT 0)) then begin
  print,'SEG_PIXVAL: Invalid column number.'
  return
endif
if ((max(row) GT s[2]-1) || (min(row) LT 0)) then begin
  print,'SEG_PIXVAL: Invalid row number.'
  return
endif

for i = 0, ncol-1 do im[col[i], row[i]] = value

MODFITS, infile, im, hdr, EXTEN_NO=segment

end

