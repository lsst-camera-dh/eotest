
pro lsst_apply_defects, infile, defectfile,  $
   DPIXPCT=dpixpct, BPIXVAL=bpixval, TRAPVAL=trapval, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose, _extra=_extra

;+
; NAME:
;   LSST_APPLY_DEFECTS
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
;    lsst_apply_defects, infile, defectfile, [, DPIXPCT=dpixpct, 
;    BPIXVAL=bpixval, TRAPVAL=trapval, PRE_COLS=pre_cols, 
;    OVER_COLS=over_cols, OVER_ROWS=over_rows, VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    DEFECTFILE: the defect mask file with the defect codes to be applied    
;
; OPTIONAL INPUTS:
;    NONE
;
; INPUT KEYWORD PARAMETERS:
;    DPIXPCT : percentage of neighborhood value for dark pixels
;    BPIXVAL : pixel value to be written for bright pixels
;    TRAPVAL : pixel value to be written for trap pixels
;    PRE_COLS:  number of prescan columns in image segment, default=10
;    OVER_COLS: number of overscan columns in image segment, default=20
;    OVER_ROWS: number of overscan rows in image segment, default=20
;    VERBOSE: if set, print some diagnostics
;    
;    'extra' keywords which may be passed to seg_apply_defects include
;    DEAD   : if set, apply 'dead' pixels
;    DARK   : if set, apply 'dark' pixels
;    BRIGHT : if set, apply 'bright' pixels
;    TRAPS  : if set, apply 'trap' pixels
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
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

dpixpct = float(dpixpct)
trappct = float(trappct)

if ((dpixpct LT 0.0) || (dpixpct GT 100.0)) then begin
  print,'LSST_APPLY_DEFECTS: Invlaid DPIXPCT. Exiting.'
  return
endif

if ((bpixval LT 0) || (bpixval GT 65535)) then begin
  print,'LSST_APPLY_DEFECTS: Invlaid BPIXVAL. Exiting.'
  return
endif

if ((trappct LT 0.0) || (trappct  GT 100.0)) then begin
  print,'LSST_APPLY_DEFECTS: Invlaid TRAPPCT. Exiting.'
  return
endif



if N_Params() NE 2 then begin
  print,'LSST_APPLY_DEFECTS: Incorrect number of parameters. Exiting.'
  return
endif

if keyword_set(VERBOSE) then verbose=1 else verbose=0

; check if the input file exists and is readable, if not exit
if (~file_test(infile, /READ)) then begin
   print,'LSST_APPLY_DEFECTS: Specified input file does not exist or is unreadable'
   return
endif

; check if the defect file exists and is readable, if not exit
if (~file_test(defectfile, /READ)) then begin
   print,'LSST_APPLY_DEFECTS: Specified defect file does not exist or is unreadable'
   return
endif

; check if the input file exists and is writable, if not exit
if (~file_test(infile, /WRITE)) then begin
   print,'LSST_APPLY_DEFECTS: Specified input file does not exist or is unwritable'
   return
endif

for segment = 1, 16 do begin
   seg_apply_defects, infile, defectfile, segment,  $
      DPIXPCT=dpixpct, BPIXVAL=bpixval, TRAPVAL=trapval, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
      VERBOSE=verbose, _extra=_extra
endfor



end

