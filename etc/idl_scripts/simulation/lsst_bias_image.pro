
pro lsst_bias_image, filename, DEFECTFILE=defectfile, $
   OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, BPIXVAL=bpixval, $ 
   SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
 ; NAME:
 ;   LSST_BIAS_CREATE
 ;
 ; PURPOSE:
 ;   Create a 16 extension FITS file simulated bias image
 ;
 ; EXPLANATION:
 ;   Create a 16 extension FITS file simulated bias image
 ;   
 ;   
 ; CALLING SEQUENCE:
 ;    LSST_BIAS_CREATE, filename [ , DEFECTFILE=defectfile, 
 ;       OFFSET=offset, GAIN=gain, RDNOISE=rdnoise, BPIXVAL=bpixval,
 ;       SEG_COLS=seg_cols, SEG_ROWS=seg_rows,
 ;       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows,
 ;       VERBOSE=verbose ]
 ;
 ; INPUTS:
 ;    filename: name of the file to create
 ;
 ; INPUT KEYWORD PARAMETERS:
 ;    DEFECTFILE : if set, use specified file to add bright defects to bias image
 ;    BPIXVAL : value to set any bright defects to
 ;    OFFSET  : image bias offset in DN
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
  print,'LSST_BIAS_IMAGE: Incorrect number of parameters. Exiting.'
  return
endif

IF N_Elements(offset)    EQ 0 THEN offset = 1000   ; image bias offset in DN
IF N_Elements(gain)      EQ 0 THEN gain = 5.0      ; system gain in e-/DN
IF N_Elements(rdnoise)   EQ 0 THEN rdnoise = 7.0   ; system read noise in e-
IF N_Elements(bpixval)   EQ 0 THEN bpixval = 65000 ; really bright pixels!
IF N_Elements(seg_cols)  EQ 0 THEN seg_cols  = 512
IF N_Elements(seg_rows)  EQ 0 THEN seg_rows  = 2002
IF N_Elements(pixval)    EQ 0 THEN pixval    = 0
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

if keyword_set(VERBOSE) then verbose = 1 else verbose = 0

gain = float(gain)
rdnoise = float(rdnoise)

lsst_image_create, filename, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PIXVAL=offset, PRE_COLS=pre_cols, OVER_COLS=over_cols, $
   OVER_ROWS=over_rows, VERBOSE=verbose
   
IF N_Elements(rdnoise) EQ 1 THEN for i = 1, 16 do $
   seg_noise, filename, i, rdnoise, GAIN=gain, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
      /full
  
IF N_Elements(defectfile) EQ 1 THEN begin
   lsst_apply_defects, filename, defectfile,  $
      BPIXVAL=bpixval,  $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
      VERBOSE=verbose, /bright
endif

end