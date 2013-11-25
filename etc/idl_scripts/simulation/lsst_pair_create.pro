
pro lsst_pair_create, filebase, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PIXVAL=pixval, IMVAL=imval, OVERVAL=overval, GAIN=gain, POISSON=poisson, $
   DEFECTFILE=defectfile, $
   NOISE=noise, PRE_COLS=pre_cols, OVER_COLS=over_cols, $
   OVER_ROWS=over_rows, VERBOSE=verbose

 ;+
 ; NAME:
 ;   LSST_PAIR_CREATE
 ;
 ; PURPOSE:
 ;   Create a pair of 16 extension FITS files such as LSST test stands do
 ;
 ; EXPLANATION:
 ;
 ; CALLING SEQUENCE:
 ;    LSST_PAIR_CREATE, filebase, [ SEG_COLS=seg_cols, SEG_ROWS=seg_rows, 
 ;       PIXVAL=pixval, IMVAL=imval, OVERVAL=overval, GAIN=gain, 
 ;       POISSON=poisson, NOISE=noise, 
 ;       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows,
 ;       CHAN_ID=chan_id, VERBOSE=verbose ]
 ;
 ; INPUTS:
 ;    filebase: base name of the files to create. Two files will be created:
 ;        filebase_flat1.fits and filebase_flat2.fits
 ;
 ; OPTIONAL INPUTS:
 ;
 ; INPUT KEYWORD PARAMETERS:
 ;
 ;    SEG_COLS: number of imaging columns in each segment, default = 512
 ;    SEG_ROWS: number of imaging rows in each segment, default = 2002
 ;    PIXVAL: value to set all pixels to unless CHAN_ID is set, default = 0
 ;    IMVAL : value to set image array pixels to, unless PIXVAL or CHAN_ID are set
 ;    GAIN: system gain in electrons per DN
 ;    POISSON : if set, add poison noise to the image area
 ;    OVERVAL: if set, the value to use for overscan pixels
 ;    NOISE: if set, the rms noise for the images
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
compile_opt IDL2

nParam = N_Params()
if nParam NE 1 then begin
  print,'LSST_PAIR_CREATE: Incorrect number of parameters. Exiting.'
  return
endif

IF N_Elements(seg_cols)  EQ 0 THEN seg_cols  = 512
IF N_Elements(seg_rows)  EQ 0 THEN seg_rows  = 2002
IF N_Elements(pixval)    EQ 0 THEN pixval    = 0
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20
IF N_Elements(noise) EQ 0 THEN noise = 0.00
IF N_Elements(gain) EQ 0 THEN gain = 1.0 ; in e-/DN

; make the first file
filename = filebase+'_flat1.fits'
print,' Creating flat image file: ',filename
lsst_image_create, filename, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
  PIXVAL=pixval, IMVAL=imval, OVERVAL=overval, GAIN=gain, POISSON=poisson, $
  DEFECTFILE=defectfile, $
  NOISE=noise, PRE_COLS=pre_cols, OVER_COLS=over_cols, $
  OVER_ROWS=over_rows, VERBOSE=verbose
  
; make the second file
filename = filebase+'_flat2.fits'
print,' Creating flat image file: ',filename
lsst_image_create, filename, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
  PIXVAL=pixval, IMVAL=imval, OVERVAL=overval, GAIN=gain, POISSON=poisson, $
  DEFECTFILE=defectfile, $
  NOISE=noise, PRE_COLS=pre_cols, OVER_COLS=over_cols, $
  OVER_ROWS=over_rows, VERBOSE=verbose

end