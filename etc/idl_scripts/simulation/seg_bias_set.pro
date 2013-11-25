
pro seg_bias_set, infile, segment, value, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, VERBOSE=verbose

;+
; NAME:
;   SEG_BIAS_SET
;
; PURPOSE:
;   Set the bias offset level in a segment to a given mean value
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    seg_bias_set, infile, segment, value [ , VERBOSE=verbose ]
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

if n_elements(PRE_COLS)  EQ 0 then pre_cols = 10
if n_elements(OVER_COLS) EQ 0 then over_cols= 20
if n_elements(OVER_ROWS) EQ 0 then over_rows= 20
if keyword_set(VERBOSE) then verbose = 1 else verbose = 0

errmsg=' '
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

if ((over_cols EQ 0) && (over_rows EQ 0)) then begin
  print,'SEG_BIAS_SET: cannot perform operation without overscan. Exiting'
  return
endif

fits_read, infile, im, hdr, exten_no=segment, /no_pdu
im = long(im) ; temporarily make image a long to do the math without clipping
s=size(im)

; find overscan offset level
if (over_cols GT 0) then bias_cols = mean(im[(s[1]-over_cols):s[1]-1, *])
if (over_rows GT 0) then bias_rows = mean(im[*,(s[2]-over_rows):s[2]-1])
if (over_cols EQ 0) then bias = long(bias_rows)
if (over_rows EQ 0) then bias = long(bias_cols)
if ((over_cols GT 0) && (over_rows GT 0)) then bias = long((bias_cols + bias_rows)/2)


; make the value a long integer
value = long(value)

diff = value - bias

; do the math
im = im + diff

; clip the result
im = im < 65535
im = im > 0

; convert back to a UINT
im = uint(im)

; save the file
MODFITS, infile, im, hdr, EXTEN_NO=segment

end

