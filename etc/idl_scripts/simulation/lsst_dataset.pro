
pro lsst_dataset, serno, OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, $
   SEG_COLS=seg_cols, SEG_ROWS=seg_rows, PRE_COLS=pre_cols, $
   OVER_COLS=over_cols, OVER_ROWS=over_rows, VERBOSE=verbose

 ;+
 ; NAME:
 ;   LSST_DATASET
 ;
 ; PURPOSE:
 ;   Create a set of simulated image data files such as LSST test stands do
 ;
 ; EXPLANATION:
 ;
 ; CALLING SEQUENCE:
 ;    LSST_DATASET, serno, [ OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, 
 ;       SEG_COLS=seg_cols, SEG_ROWS=seg_rows, PRE_COLS=pre_cols, 
 ;       OVER_COLS=over_cols, OVER_ROWS=over_rows, VERBOSE=verbose ]
 ;
 ; INPUTS:
 ;    serno: string representing the serial number of the device e.g. '000_00'
 ;
 ; OPTIONAL INPUTS:
 ;
 ; INPUT KEYWORD PARAMETERS:
 ;
 ;    OFFSET  : image bias offset in DN
 ;    GAIN    : system gain, e-/DN
 ;    RDNOISE : system read noise in e-
 ;    SEG_COLS: number of imaging columns in each segment, default = 512
 ;    SEG_ROWS: number of imaging rows in each segment, default = 2002
 ;    PRE_COLS:  number of prescan columns in image segment, default=10
 ;    OVER_COLS: number of overscan columns in image segment, default=20
 ;    OVER_ROWS: number of overscan rows in image segment, default=20
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

if N_Params() NE 1 then begin
  print,'LSST_DATASET: Incorrect number of parameters. Exiting.'
  return
endif

IF N_Elements(offset)    EQ 0 THEN offset = 1000   ; image bias offset in DN
IF N_Elements(gain)      EQ 0 THEN gain = 5.0      ; system gain in e-/DN
IF N_Elements(rdnoise)   EQ 0 THEN rdnoise = 7.0 ; system read noise in e-
IF N_Elements(seg_cols)  EQ 0 THEN seg_cols  = 512
IF N_Elements(seg_rows)  EQ 0 THEN seg_rows  = 2002
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20
IF N_Elements(gain) EQ 0 THEN gain = 1.0 ; in e-/DN

if keyword_set(VERBOSE) then verbose=1 else verbose=0

; create a mask image for use in the set
maskfile=string(serno,'_mask.fits',format='(A,A)')
lsst_mask_image, maskfile, EDGE=9, MIDLINE=5, $
   SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose
   
; create a simulated defect file 
defectfile=string(serno,'_defect.fits',format='(A,A)')
lsst_defect_image, defectfile, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

; create a results file for use in analysis
resultsfile=string(serno,'_result.fits',format='(A,A)')
lsst_image_create, resultsfile, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PIXVAL=0, PRE_COLS=pre_cols, OVER_COLS=over_cols, $
   OVER_ROWS=over_rows, VERBOSE=verbose

; create dark current images and associated biases.
for i = 1, 5 do begin
   ; make bias image
   filename = string(serno,'_dark_bias_',i,'.fits',format='(A,A,I03,A)')
   lsst_bias_image, filename, DEFECTFILE=defectfile, $
      OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, BPIXVAL=bpixval, $ 
      SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
      
   ; then a dark, which is a flat at very low "flux" for 500 seconds
   ; darks don't get 'dark' defects, so do bright defects after making the image
   filename = string(serno,'_dark_dark_',i,'.fits',format='(A,A,I03,A)')
   lsst_flat_image, filename, 0.1, 500, $
      OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, $
      SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
   lsst_apply_defects, filename, defectfile, BPIXVAL=65000, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
      VERBOSE=verbose, /bright, /dead, /traps
endfor

; create a big stack o' pairs of flats with different exposure times
flux = 10000  ; e-/pixel/second
exptimes = [ 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00, 3.00, $
             4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 11.0, 12.0, $
             13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0 ]
print,exptimes
for i = 1, 25 do begin
  filebase=string(serno,'_flat_',i,format='(A,A,I03)')
  lsst_flat_pair, filebase, flux, exptimes[i], $
     DEFECTFILE=defectfile, BPIXVAL=bpixval, $
     OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, $
     SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
     PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
     VERBOSE=verbose
endfor 
               
; make some fe55 x-ray images
for i = 1, 5 do begin
   ; make bias image
   filename = string(serno,'_fe55_',i,'.fits',format='(A,A,I03,A)')
   lsst_bias_image, filename, DEFECTFILE=defectfile, $
      OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, BPIXVAL=bpixval, $ 
      SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows
             
   for seg = 1, 16 do begin
      seg_fe55_add, filename, seg, GAIN=gain, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
      VERBOSE=verbose
   endfor
   
endfor            
end