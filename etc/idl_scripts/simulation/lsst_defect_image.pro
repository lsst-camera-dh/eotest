

pro lsst_defect_image, filename, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
 ; NAME:
 ;   LSST_DEFECT_IMAGE
 ;
 ; PURPOSE:
 ;   Create a 16 extension FITS file defect mask image with defect pixels 
 ;   identified via bit mask. If no bits are set the pixel is neither masked
 ;   nor a defect. Bits are defined in the following manner:
 ;   
 ;   0 : masked pixel, masked for any reason, not used in calculations
 ;   1 : dead pixel: pixel does not respond to light
 ;   2 : dark pixel: pixel has less response to light vs other pixels
 ;   3 : dark column : pixel in a column defined as dark due to dark defects
 ;   4 : bright pixel : pixel has excess dark current
 ;   5 : bright column : pixel in a column defined as bright due to bright defects
 ;   6 : trap pixel : a pixel that traps charge in excess of some limit
 ;   7 : trap column : pixel in a column that contains a trap
 ;   etc. 
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
 ; EXPLANATION:
 ;
 ; CALLING SEQUENCE:
 ;    LSST_DEFECT_IMAGE, filename, [ SEG_COLS=seg_cols, SEG_ROWS=seg_rows,
 ;       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows,
 ;       VERBOSE=verbose ]
 ;
 ; INPUTS:
 ;    filename: name of the file to create
 ;
 ; OPTIONAL INPUTS:
 ;
 ; INPUT KEYWORD PARAMETERS:
 ;
 ;    SEG_COLS: number of imaging columns in each segment, default = 512
 ;    SEG_ROWS: number of imaging rows in each segment, default = 2002
 ;    PRE_COLS:  number of prescan columns in image segment, default=10
 ;    OVER_COLS: number of overscan columns in image segment, default=20
 ;    OVER_ROWS: number of overscan rows in image segment, default=20
 ;    CHAN_ID: if set, make all pixels in extension equal to extension number, not pixval
 ;    VERBOSE: if set, print some diagnostics
 ;    
 ; OUTPUT KEYWORD PARAMETERS:
 ;    NONE
 ;
 ; EXAMPLES:
 ;
 ;
 ; HISTORY:
 ; Written by: P. Doherty  April, 2013
 ;
 ;-
 ;-----------------------------------------------------------------------------   
compile_opt IDL2

nParam = N_Params()
if nParam NE 1 then begin
  print,'LSST_DEFECT_IMAGE: Incorrect number of parameters. Exiting.'
  return
endif

IF N_Elements(seg_cols)  EQ 0 THEN seg_cols  = 512
IF N_Elements(seg_rows)  EQ 0 THEN seg_rows  = 2002
IF N_Elements(pixval)    EQ 0 THEN pixval    = 0
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

if keyword_set(VERBOSE) then verbose = 1 else verbose = 0

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

if (verbose) then print,'Making mask image...'
lsst_mask_image, filename, EDGE=5, MIDLINE=5, $
   SEG_COLS=seg_cols, SEG_ROWS=seg_rows, PRE_COLS=pre_cols, OVER_COLS=over_cols, $
   OVER_ROWS=over_rows,  VERBOSE=verbose

; how about 50 dark pixels in each segment?
x = uint(indgen(50) * (seg_cols-20)/50.0 + pre_cols)
y = uint(indgen(50) * (seg_rows-20)/50.0)
if (verbose) then print,'Adding 50 dark pixels...'
for segment = 1, 16 do $
  seg_add_defect, filename, segment, x, y, 'dpix', $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows

; and 50 bright pixels in each segment?
if (verbose) then print,'Adding 50 bright pixels...'
x = seg_cols-x-2
for segment = 1, 16 do $
   seg_add_defect, filename, segment, x, y, 'bpix', $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
      
; how about 1 full bright column defect and two partials (one less than 20)?
if (verbose) then print,'Adding 3 bright columns...'
for segment = 1, 16 do $
   seg_colval, filename, segment, 105, bpix, MIN_ROW=0, MAX_ROW=seg_rows, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
for segment = 1, 16 do $
   seg_colval, filename, segment, 220, bpix, MIN_ROW=400, MAX_ROW=800,   $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
for segment = 1, 16 do $
   seg_colval, filename, segment, 450, bpix, MIN_ROW=1000, MAX_ROW=1015, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
  
; how about 1 full dark column defect and 2 partials, one less than 100?
if (verbose) then print,'Adding 3 dark columns...'
for segment = 1, 16 do $
   seg_colval, filename, segment, 115, dpix, MIN_ROW=0, MAX_ROW=seg_rows, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
for segment = 1, 16 do $
   seg_colval, filename, segment, 230, dpix, MIN_ROW=400, MAX_ROW=800,   $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
for segment = 1, 16 do $
   seg_colval, filename, segment, 460, dpix, MIN_ROW=1000, MAX_ROW=1030, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows

; a random smattering of dead pixels?
if (verbose) then print,'Adding dead pixels...'
for segment = 1, 16 do begin
     col = uint(randomu(seed,20)*seg_cols) + pre_cols
     row = uint(randomu(seed,20)*seg_rows)
     seg_add_defect, filename, segment, col, row, 'dead', $
         PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
endfor

end