
pro lsst_fe55_image, filename, DEFECTFILE=defectfile, BPIXVAL=bpixval, $
   OFFSET=offset, SIGNAL=signal, GAIN=gain, RDNOISE=rdnoise, POISSON=poisson, $
   SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
 ; NAME:
 ;   LSST_FLAT_IMAGE
 ;
 ; PURPOSE:
 ;   Create a 16 extension FITS file simulated flat field image
 ;
 ; EXPLANATION:
 ;   Create a 16 extension FITS file simulated flat field image
 ;   
 ;   
 ; CALLING SEQUENCE:
 ;    LSST_FLAT_IMAGE, filename [ , OFFSET=offset, GAIN=gain, 
 ;       RDNOISE=rdnoise, POISSON=poisson, 
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
 ;    DEFECTFILE : if set, use specified file to add bright defects to bias image
 ;    BPIXVAL : value to set any bright defects to
 ;    OFFSET  : image bias offset in DN
 ;    SIGNAL  : flat field signal level, in electrons
 ;    GAIN    : system gain, e-/DN
 ;    RDNOISE : system read noise in e-
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
  print,'LSST_FLAT_IMAGE: Incorrect number of parameters. Exiting.'
  return
endif

IF N_Elements(offset)    EQ 0 THEN offset = 1000   ; image bias offset in DN
IF N_Elements(signal)    EQ 0 THEN signal = 0      ; flat field signal, e-
IF N_Elements(gain)      EQ 0 THEN gain = 5.0      ; system gain in e-/DN
IF N_Elements(rdnoise)   EQ 0 THEN readnoise = 7.0 ; system read noise in e-
IF N_Elements(bpixval)   EQ 0 THEN bpixval = 65000 ; really bright pixels!
IF N_Elements(seg_cols)  EQ 0 THEN seg_cols  = 512
IF N_Elements(seg_rows)  EQ 0 THEN seg_rows  = 2002
IF N_Elements(pixval)    EQ 0 THEN pixval    = 0
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20
if keyword_set(POISSON) then poisson = 1 else poisson = 0
if keyword_set(VERBOSE) then verbose = 1 else verbose = 0

gain = float(gain)
rdnoise = float(rdnoise)

;lsst_image_create, filename, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
;   PIXVAL=offset, $
;   PRE_COLS=pre_cols, OVER_COLS=over_cols, $
;   OVER_ROWS=over_rows, VERBOSE=verbose

lsst_bias_image, filename, $
   OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, $ 
   SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

; add simultaed fe55 events here
; 
;for i = 1, 16 do begin
;   seg_add_signal, filename, i, 5000, GAIN=gain, POISSON=poisson, $
;      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
;      VERBOSE=verbose
;endfor
  
IF N_Elements(defectfile) EQ 1 THEN begin
lsst_apply_defects, filename, defectfile,  $
   BPIXVAL=65000,  $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose, /bright
endif



end