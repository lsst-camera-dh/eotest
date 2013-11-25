
pro seg_add_defect, infile, segment, col, row, defect, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
; NAME:
;   SEG_ADD_DEFECT
;
; PURPOSE:
;   Put a defect code in place of a single pixel in a FITS file extension. 
;
; EXPLANATION:
;   defect pixels are identified via bit mask. If no bits are set the pixel is 
;   neither masked nor a defect. Bits are defined in the following manner:
;   
;   0 : masked pixel, masked for any reason, not used in calculations
;   1 : dead pixel: pixel does not respond to light
;   2 : dark pixel: pixel has less response to light vs other pixels
;   3 : dark column : pixel in a column defined as dark due to dark defects
;   4 : bright pixel : pixel has excess dark current
;   5 : bright column : pixel in a column defined as bright due to bright defects
;   6 : trap pixel : a pixel that traps charge in excess of some limit
;   7 : trap column : pixel in a column that contains a trap
;  
;   Pixels may fall into multiple categories. Some defect examples:
;   
;   Masked pixel : 0x1
;   Dead Pixel : 0x2
;   Dark pixel not in a dark column : 0x4
;   Functional pixel in a dark column : 0x8
;   Dark pixel in a dark column : 0x4 || 0x8 = 0xc
;   etc.
;   
;   The type of defect, 'defect', added to the image may be set either as a 
;   hexadecimal number or as a string from the following list: 
;   
;   'mask' = '01'XU = masked pixel
;   'dead' = '02'XU = dead pixel
;   'dpix' = '04'XU = dark pixel
;   'dcol' = '08'XU = dark column
;   'bpix' = '10'XU = bright pixel
;   'bcol' = '20'XU = bright column
;   'trap' = '40'XU = trap pixel
;   'tcol' = '80'XU = trap column
;
;
; CALLING SEQUENCE:
;    SEG_ADD_DEFECT, infile, segment, col, row, defect [ , PRE_COLS=pre_cols, 
;        OVER_COLS=over_cols, OVER_ROWS=over_rows, VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;    COL: column number of the pixel to set, may be an array
;    ROW: row number of the pixel to set, may be an array
;    DEFECT: defect code or defect type string
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
;    NONE
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

if keyword_set(VERBOSE) then verbose=1 else verbose=0

if N_Params() NE 5 then begin
  print,'SEG_ADD_DEFECT: Incorrect number of parameters. Exiting.'
  return
endif

;if (n_elements(col) NE n_elements(row)) then begin
;  print,'SEG_ADD_DEFECT: Mismatched row and column arrays. Exiting.'
;  return
;endif
  
ncol = n_elements(col)
nrow = n_elements(row)

errmsg='No error.'
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

;   0 : masked pixel, masked for any reason, not used in calculations
;   1 : dead pixel: pixel does not respond to light
;   2 : dark pixel: pixel has less response to light vs other pixels
;   3 : dark column : pixel in a column defined as dark due to dark defects
;   4 : bright pixel : pixel has excess dark current
;   5 : bright column : pixel in a column defined as bright due to bright defects
;   6 : trap pixel : a pixel that traps charge in excess of some limit
;   7 : trap column : pixel in a column that contains a trap

mask = '01'XU
dead = '02'XU
dpix = '04'XU
dcol = '08'XU
bpix = '10'XU
bcol = '20'XU
trap = '40'XU
tcol = '80'XU

; check the defect parameter type
deftype = size(defect,/tn)

case deftype of
   'INT'    : value = uint(defect)  ; if an integer, just use that
   'UINT'   : value = defect   ; if an unsigned integer, just use that
   'STRING' : begin
       case defect of
          'mask' : begin if (verbose) then print, 'mask pixel'    & value = mask & end
          'mcol' : begin if (verbose) then print, 'mask column'   & value = mask & end
          'mrow' : begin if (verbose) then print, 'mask row'      & value = mask & end
          'dead' : begin if (verbose) then print, 'dead pixel'    & value = dead & end
          'dpix' : begin if (verbose) then print, 'dark pixel'    & value = dpix & end
          'dcol' : begin if (verbose) then print, 'dark column'   & value = dcol & end
          'bpix' : begin if (verbose) then print, 'bright pixel'  & value = bpix & end
          'bcol' : begin if (verbose) then print, 'bright column' & value = bcol & end
          'trap' : begin if (verbose) then print, 'trap pixel'    & value = trap & end
          'tcol' : begin if (verbose) then print, 'trap column'   & value = tcol & end
          'mask pixel'    : begin if (verbose) then print, 'mask pixel'    & value = mask & end
          'mask column'   : begin if (verbose) then print, 'mask column'   & value = mask & end
          'mask row'      : begin if (verbose) then print, 'mask row'      & value = mask & end
          'dead pixel'    : begin if (verbose) then print, 'dead pixel'    & value = dead & end
          'dark pixel'    : begin if (verbose) then print, 'dark pixel'    & value = dpix & end
          'dark column'   : begin if (verbose) then print, 'dark column'   & value = dcol & end
          'bright pixel'  : begin if (verbose) then print, 'bright pixel'  & value = bpix & end
          'bright column' : begin if (verbose) then print, 'bright column' & value = bcol & end
          'trap pixel'    : begin if (verbose) then print, 'trap pixel'    & value = trap & end
          'trap column'   : begin if (verbose) then print, 'trap column'   & value = tcol & end
       else: begin
          print,'SEG_ADD_DEFECT: Invalid defect type string'
          return
          end
       endcase
    end
else: begin
   print,'SEG_ADD_DEFECT: Invalid defect type'
   return
   end
endcase

; read in the image segment
fits_read, infile, im, hdr, exten_no=segment
s = size(im)

; check for valid column numbers, defects only exist in image array, but overscan may be masked
if (value EQ mask) then begin
   if ((max(col) GT s[1]-1) || (min(col) LT 0)) then begin
      print,'SEG_ADD_DEFECT: Invalid column number.'
      return
   endif
endif else begin
   if ((max(col) GT s[1]-over_cols) || (min(col) LT pre_cols)) then begin
      print,'SEG_ADD_DEFECT: Invalid column number.'
      return
   endif
endelse

; check for valid row number, defects only exist in image array, but overscan may be masked
 if (value EQ mask) then begin
   if ((max(row) GT s[2]-1) || (min(row) LT 0)) then begin
      print,'SEG_ADD_DEFECT: Invalid row number.'
      return
   endif
endif else begin
   if ((max(row) GT s[2]-over_rows) || (min(row) LT 0)) then begin
      print,'SEG_ADD_DEFECT: Invalid row number.'
      return
   endif
endelse

; if user specified an integer value, do what they asked to single pixels
if ((deftype EQ 'INT') || (deftype EQ 'UINT')) then begin
  if (n_elements(col) NE n_elements(row)) then begin
     print,'SEG_ADD_DEFECT: Mismatched row and column arrays. Exiting.'
     return
  endif
  for i = 0, ncol-1 do im[col[i],row[i]] = im[col[i],row[i]] OR value
  if file_test(infile, /WRITE) then MODFITS, infile, im, hdr, EXTEN_NO=segment
  return
endif

; if mask column, mask whole columns including overscan rows, ignore row parameter
if (defect EQ 'mcol') then $
   for i = 0, ncol-1 do im[col[i],*] = im[col[i],*] OR value

; if dark col, bright col, or trap col, do whole column in image array only
if (value EQ dcol) || (value EQ bcol) || (value EQ tcol) then $
   for i = 0, ncol-1 do im[col[i], 0:s[2]-over_rows] = im[col[i], 0:s[2]-over_rows] OR value

; if mask row, mask the whole row
if (defect EQ 'mrow') then $
   for i = 0, nrow-1 do im[*,row[i]] = im[*,row[i]] OR value

; for single pixel mask or defect, just do one pixel
if ((defect EQ 'mask') || (defect EQ 'mask pixel') || (value EQ dpix) || $
    (value EQ bpix) || (defect EQ 'trap') || (defect EQ 'trap pixel') || $
    (value EQ dead)) then begin 
      if (n_elements(col) NE n_elements(row)) then begin
         print,'SEG_ADD_DEFECT: Mismatched row and column arrays. Exiting.'
         return
      endif
      for i = 0, ncol-1 do $
         im[col[i], row[i]] = im[col[i], row[i]] OR value 
endif
; save the result
if file_test(infile, /WRITE) then MODFITS, infile, im, hdr, EXTEN_NO=segment

end

