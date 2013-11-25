
pro seg_const_add, infile, segment, value, VERBOSE=verbose

;+
; NAME:
;   SEG_CONST_ADD
;
; PURPOSE:
;   Add a constant value to every pixel in an image segment
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    seg_const_add, infile, segment, value [ , VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    VALUE: value to add to the pixels
;
; OPTIONAL INPUTS:
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

errmsg=' '
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

fits_read, infile, im, hdr, exten_no=segment
im = long(im) ; temporarily make image a long to do the math without clipping

; make the value a long integer
value = long(value)

; do the math
im = im + value

; clip the result
im = im < 65535
im = im > 0

; convert back to a UINT
im = uint(im)

; save the file
MODFITS, infile, im, hdr, EXTEN_NO=segment

end

