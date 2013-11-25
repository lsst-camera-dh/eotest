
pro seg_apply_defects, infile, defectfile, segment,  $
   DEAD=dead, DARK=DARK, BRIGHT=bright, TRAPS=traps, $
   DPIXPCT=dpixpct, BPIXVAL=bpixval, TRAPVAL=trapval, BOX_SIZE=box_size, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
; NAME:
;   SEG_APPLY_DEFECTS
;
; PURPOSE:
;   Add defects to an image segment using a defect mask image as guide 
;
; EXPLANATION:
;   In the defect mask image defect pixels are identified via bit mask. 
;   If no bits are set the pixel is neither masked nor a defect. Bits are 
;   defined in the following manner:
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
;   Pixels identified as defects in the defect file are set to appropriate values 
;   in the input file. The values to use for each type of defect can be set via
;   input keywords or the default values my be used. 
;   
;   'Dead' pixels are replaced with the mean value found in the overscan. 
;   'Dark' pixels are made to be some percentage less bright than their neighbors.
;   'Bright' pixels are set to a defined constant value
;   'Trap' pixels are made to be some percentage less bright than the previous
;       pixel in the column.
;   
; CALLING SEQUENCE:
;    seg_apply_defects, infile, defectfile, segment, [, DEAD=dead, 
;      DEAD=dead, DARK=DARK, BRIGHT=bright, TRAPS=traps,
;      BPIXVAL=bpixval, TRAPVAL=trapval, BOX_SIZE=box_size, PRE_COLS=pre_cols, 
;      OVER_COLS=over_cols, OVER_ROWS=over_rows, VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    DEFECTFILE: the defect mask file with the defect codes to be applied    
;    SEGMENT: number of the extension to operate on
;
; OPTIONAL INPUTS:
;    NONE
;
; INPUT KEYWORD PARAMETERS:
;    DEAD   : if set, apply 'dead' pixels
;    DARK   : if set, apply 'dark' pixels
;    BRIGHT : if set, apply 'bright' pixels
;    TRAPS  : if set, apply 'trap' pixels
;    DPIXPCT : percentage of neighborhood value for dark pixels
;    BPIXVAL : pixel value to be written for bright pixels
;    TRAPVAL : pixel value to be written for trap pixels
;    BOX_SIZE: size of the box to use when computing dark pixel value
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

IF N_Elements(dpixpct)  EQ 0 THEN dpixpct  = 50.0
IF N_Elements(bpixval)  EQ 0 THEN bpixval  = 65000
IF N_Elements(trappct)  EQ 0 THEN trappct  = 50.0
IF N_Elements(box_size)  EQ 0 THEN box_size=5

dpixpct = float(dpixpct)
trappct = float(trappct)

if ((dpixpct LT 0.0) || (dpixpct GT 100.0)) then begin
  print,'SEG_APPLY_DEFECTS: Invlaid DPIXPCT. Exiting.'
  return
endif

if ((bpixval LT 0) || (bpixval GT 65535)) then begin
  print,'SEG_APPLY_DEFECTS: Invlaid BPIXVAL. Exiting.'
  return
endif

if ((trappct LT 0.0) || (trappct  GT 100.0)) then begin
  print,'SEG_APPLY_DEFECTS: Invlaid TRAPPCT. Exiting.'
  return
endif

IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

if N_Params() NE 3 then begin
  print,'SEG_APPLY_DEFECTS: Incorrect number of parameters. Exiting.'
  return
endif

if keyword_set(VERBOSE) then verbose=1 else verbose=0

errmsg='No error.'
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

; read in the input file and compute boundaries
fits_read, infile, im, hdr, exten_no=segment
s=size(im)
min_col = pre_cols              ; stay out of prescan
min_row = 0                     ; bottom row is fine
max_col = s[1] - (over_cols+1)  ; stay out of the overscan in X
max_row = s[2] - (over_rows+1)  ; stay out of the overscan in Y

; read in the input file and compute boundaries
fits_read, defectfile, defects, dhdr, exten_no=segment
s2=size(defects)

if ((s[1] NE s2[1]) || (s[2] NE s2[2])) then begin
  print,'SEG_APPLY_DEFECTS: Input file sizes do not match. Exiting.'
  return
endif
  
; bits used to define defects in the defect mask file
;   0 : masked pixel, masked for any reason, not used in calculations
;   1 : dead pixel: pixel does not respond to light
;   2 : dark pixel: pixel has less response to light vs other pixels
;   3 : dark column : pixel in a column defined as dark due to dark defects
;   4 : bright pixel : pixel has excess dark current
;   5 : bright column : pixel in a column defined as bright due to bright defects
;   6 : trap pixel : a pixel that traps charge in excess of some limit
;   7 : trap column : pixel in a column that contains a trap

mask_code = '01'XU
dead_code = '02'XU
dpix_code = '04'XU
dcol_code = '08'XU
bpix_code = '10'XU
bcol_code = '20'XU
trap_code = '40'XU
tcol_code = '80'XU


; set dead pixels to bias level found in overscan
if keyword_set(DEAD) then begin
   bias = uint(mean(im[s[1]-over_cols:s[1]-1, 0:s[2]-1]))
   im[where(defects EQ dead_code)] = uint(bias)
endif

; set dark pixels to a specified fraction of their nearest neighbors
if keyword_set(DARK) then begin
   darks = where(defects EQ dpix_code, count)
   if count GT 0 then begin
      x = darks MOD s[1]
      y = darks / s[1]
      for i = 0, count-1 do begin
         ; buld a box that doesn't go over the edges of the imaging area
         xmin = (x[i] - box_size/2) > uint(min_col)   ; can't be outside the image negative or in prescan
         ymin = (y[i] - box_size/2) > uint(min_row)   ; can't be outside the image negative
         xmax = (xmin + box_size) < max_col        ; can't go into overscan in X
         ymax = (ymin + box_size) < max_row        ; can't go into overscan in Y
         xmin = xmax - box_size                    ; that's where it has to be in X
         ymin = ymax - box_size                    ; that's where it has to be in Y
         ; extract relevant sections from image and mask
         im_box = im[xmin:xmax,ymin:ymax]
         defect_box = defects[xmin:xmax,ymin:ymax]
         ; calculate median value in im_box excluding any defects marked in defect_box
         box_median = median(im_box[where(defect_box EQ 0)]) 
         ; set the pixel value to dpixpct of the surrounding median value
         im[x[i],y[i]] = box_median * (dpixpct/100.0)
      endfor
   endif
endif

; set bright pixels to the specified value
if keyword_set(BRIGHT) then begin
   im[where(defects EQ bpix_code)] = uint(bpixval)
endif

; set traps to something
if keyword_set(TRAPS) then begin
;   trappix = where(defects EQ trap_code, count)
;   if count GT 0 then begin
;      x = trappix MOD s[1]
;      y = trappix / s[1]
;      for i = 0, count-1 do begin
;         ; buld a box that doesn't go over the edges of the imaging area
;         xmin = (x[i] - box_size/2) > uint(min_col)   ; can't be outside the image negative or in prescan
;         ymin = (y[i] - box_size/2) > uint(min_row)   ; can't be outside the image negative
;         xmax = (xmin + box_size) < max_col        ; can't go into overscan in X
;         ymax = (ymin + box_size) < max_row        ; can't go into overscan in Y
;         xmin = xmax - box_size                    ; that's where it has to be in X
;         ymin = ymax - box_size                    ; that's where it has to be in Y
;         ; extract relevant sections from image and mask
;         im_box = im[xmin:xmax,ymin:ymax]
;         defect_box = defects[xmin:xmax,ymin:ymax]
;         ; calculate median value in im_box excluding any defects marked in defect_box
;         box_median = median(im_box[where(defect_box EQ 0)]) 
;         ; set the pixel value to dpixpct of the surrounding median value
;         im[x[i],y[i]] = box_median * (trappct/100.0)
;      endfor
;   endif
endif

; save the result
MODFITS, infile, im, hdr, EXTEN_NO=segment

end

