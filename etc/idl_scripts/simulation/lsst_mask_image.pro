
pro lsst_mask_image, filename, EDGE=edge, MIDLINE=midline, $
   SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
 ; NAME:
 ;   LSST_MASK_IMAGE
 ;
 ; PURPOSE:
 ;   Create a 16 extension FITS file defect mask image
 ;
 ; EXPLANATION:
 ;   Create a 16 extension FITS file mask image with masked pixels 
 ;   identified via bit mask. If bit 0 is set, the bit is masked and will
 ;   not be used in defect identification by later routines that use the 
 ;   mask image. If no bits are set the pixel is not masked. All pixels 
 ;   are initially set to zero. Overscan pixels are then set to 1 (masked).
 ;   Optionally, the edges of the device and the area near the midline blooming
 ;   stop may also be masked.
 ;   
 ; CALLING SEQUENCE:
 ;    LSST_MASK_IMAGE, filename [ , EDGE=edge, MIDLINE=midline, 
 ;       SEG_COLS=seg_cols, SEG_ROWS=seg_rows,
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
 ;    EDGE: number of pixels around edge of device to mask off. Default = 0
 ;    MIDLINE: number of pixels at midline of device to mask off. Default = 0
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

if N_Params() NE 1 then begin
  print,'LSST_MASK_IMAGE: Incorrect number of parameters. Exiting.'
  return
endif

IF N_Elements(edge)      EQ 0 THEN edge = 0         ; number of edge columns to mask
IF N_Elements(midline)   EQ 0 THEN midline  = 0     ; number of midline rows to mask
IF N_Elements(seg_cols)  EQ 0 THEN seg_cols  = 512  
IF N_Elements(seg_rows)  EQ 0 THEN seg_rows  = 2002
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

if keyword_set(VERBOSE) then verbose = 1 else verbose = 0

if (verbose) then print,'Making basic image...'
lsst_image_create, filename, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PIXVAL=0, PRE_COLS=pre_cols, OVER_COLS=over_cols, $
   OVER_ROWS=over_rows,  VERBOSE=verbose

;; mask off the overscan areas, we don't look for defects there
if (verbose) then print,'Masking overscan...'
for segment = 1, 16 do begin
    seg_overval, filename, segment, 1, $
       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
       VERBOSE=verbose
endfor

; optionally mask off the edges of the array, outside columns
; of segments 1, 8, 9, and 16, and the bottom rows of each segment
if (edge NE 0) then begin
   if (verbose) then print,'Masking edges...'
   cols=indgen(edge) + (pre_cols + seg_cols - edge)
   seg_add_defect, filename, 1, cols, 0, 'mcol', $
       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows

   cols = indgen(edge) + pre_cols
   seg_add_defect, filename, 8, cols, 0, 'mcol', $
       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows

   cols = indgen(edge) + (pre_cols + seg_cols - edge)
   seg_add_defect, filename, 9, cols, 0, 'mcol', $
       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows

   cols = indgen(edge) + pre_cols
   seg_add_defect, filename, 16, cols, 0, 'mcol', $
       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows

   for segment = 1, 16 do begin
      rows=indgen(edge)
      seg_add_defect, filename, segment, 0, rows, 'mrow', $
          PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
   endfor
endif

; optionally include a mask of the rows on either side of the midline
; blooming stop in all segments
if (midline NE 0) then begin
   if (verbose) then print,'Masking midline...'
   ; mask off the last imaging rows of each segment (blooming stop)
   for segment = 1, 16 do begin
      rows=indgen(midline) + seg_rows - midline
      seg_add_defect, filename, segment, 0, rows, 'mrow', $
          PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
   endfor
endif

end