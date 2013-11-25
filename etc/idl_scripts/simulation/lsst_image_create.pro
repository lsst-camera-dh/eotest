
pro lsst_image_create, filename, SEG_COLS=seg_cols, SEG_ROWS=seg_rows, $
   PIXVAL=pixval, IMVAL=imval, OVERVAL=overval, DATATYPE=datatype, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, $
   OVER_ROWS=over_rows, CHAN_ID=chan_id, VERBOSE=verbose

 ;+
 ; NAME:
 ;   LSST_IMAGE_CREATE
 ;
 ; PURPOSE:
 ;   Create a 16 extension FITS file such as LSST test stands do
 ;
 ; EXPLANATION:
 ;
 ; CALLING SEQUENCE:
 ;    LSST_IMAGE_CREATE, filename, [ SEG_COLS=seg_cols, SEG_ROWS=seg_rows, 
 ;       PIXVAL=pixval, IMVAL=imval, OVERVAL=overval,
 ;       PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, 
 ;       CHAN_ID=chan_id, VERBOSE=verbose ]
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
 ;    PIXVAL: value to set all pixels to unless CHAN_ID is set, default = 0
 ;    IMVAL : value to set image array pixels to, unless PIXVAL or CHAN_ID are set
 ;    OVERVAL: if set, the value to use for overscan pixels, unless PIXVAL or CHAN_ID are set
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
  print,'LSST_IMAGE_CREATE: Incorrect number of parameters. Exiting.'
  return
endif

IF N_Elements(seg_cols)  EQ 0 THEN seg_cols  = 512
IF N_Elements(seg_rows)  EQ 0 THEN seg_rows  = 2002
IF N_Elements(pixval)    EQ 0 THEN pixval    = 0
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20
IF N_Elements(datatype)  EQ 0 THEN datatype  = 'uint'
if keyword_set(VERBOSE) then verbose=1 else verbose=0


CASE datatype OF
   'byte'  : dtype=1
   'uint'  : dtype=2
   'long'  : dtype=3
   'float' : dtype=4
   'double': dtype=5
   ELSE: begin
    print, "LSST_IMAGE_CREATE: undefined data type. Use 'byte', 'uint', 'long', or 'float'."
    print, 'Exiting."
    end
ENDCASE

; open the output file and make a simple primary header of the correct type
FITS_OPEN, filename, fcb, /write
mkhdr, phdr, dtype, /extend

detcols = 8 * (seg_cols + pre_cols + over_cols)
detrows = 2 * (seg_rows + over_rows)
if keyword_set(VERBOSE) then print,'detrows = ',detrows, '   detcols = ', detcols
SXADDPAR, phdr, 'DET_COLS', detcols,  'Number of imaging columns per image',   FORMAT='(I)'
SXADDPAR, phdr, 'DET_ROWS', detrows,  'Number of imaging rows per image',      FORMAT='(I)'
SXADDPAR, phdr, 'SEG_COLS', seg_cols, 'Number of imaging columns per segment', FORMAT='(I)'
SXADDPAR, phdr, 'SEG_ROWS', seg_rows, 'Number of imaging rows per segment',    FORMAT='(I)'
SXADDPAR, phdr, 'PRE_COLS', pre_cols, 'Number of prescan columns per segment', FORMAT='(I)'
SXADDPAR, phdr, 'OVR_COLS', over_cols,'Number of overscan columns per segment',FORMAT='(I)'
SXADDPAR, phdr, 'OVR_ROWS', over_rows,'Number of overscan rows per segment',   FORMAT='(I)'
SXADDPAR, phdr, 'ORIGIN  ', 'Simulation', 'Source of the image',               FORMAT='(A)'
SXADDPAR, phdr, 'INSTRUME', 'None', 'Instrument used to colect data',          FORMAT='(A)'

detstr=string('[0:',detcols,',0:',detrows,']',format='(A,I0,A,I0,A)')
SXADDPAR, phdr, 'DETSIZE ', detstr, FORMAT='(A)'
SXADDPAR, phdr, 'CCD_TYPE', 'LSST',  'CCD type', FORMAT='(A)'
SXADDPAR, phdr, 'CCD_SER ', '000-00','CCD serial number', FORMAT='(A)'
SXADDPAR, phdr, 'BINNING ', '1x1',   '[pixelX x pixelY] chip binning', FORMAT='(A)'
SXADDPAR, phdr, 'BINX    ',  1 ,     '[pixels] binning along X axis', FORMAT='(I)'
SXADDPAR, phdr, 'BINY    ',  1 ,     '[pixels] binning along X axis', FORMAT='(I)'

; write the primary FITS header
FITS_WRITE, fcb, 0, phdr

;
; primary header complete, fill in 16 headers and data arrays

; make a data array to store in each of the extensions
segw = uint(seg_cols + pre_cols + over_cols)
segh = uint(seg_rows + over_rows)
im = uintarr(segw, segh)
im[*] = uint(pixval) ; set all pixels to pixval, may update at end

; create a header for the segment image
MKHDR, seghdr, im, /IMAGE 

; add keywords that are common to all segments
SXADDPAR, seghdr, 'BZERO   ', 32768, FORMAT='(I)', AFTER='GCOUNT'
SXADDPAR, seghdr, 'BSCALE  ', 1, FORMAT='(I)', AFTER='BZERO'
SXADDPAR, seghdr, 'O_BZERO ', 32768, 'Original Data is Unsigned Integer', FORMAT='(I)', AFTER='BZERO'

; create and add each of the segments in turn
for i = 1, 16 do begin
   SXADDPAR, seghdr, 'CHANNEL ', i, 'channel number', FORMAT='(I)', AFTER='O_BZERO'
   SXADDPAR, seghdr, 'DETSIZE ', detstr, FORMAT='(A)', AFTER='CHANNEL'
   datasec=string('[1:',segw,',1:',segh,']',format='(A,I0,A,I0,A)')
   SXADDPAR, seghdr, 'DATASEC ', datasec, FORMAT='(A)',AFTER='DETSIZE'

   ; do DETSEC correctly here, this *seems to be* correct.
   if (i LT 9) then begin
       x1 = (i-1) * segw
       x2 = x1 + segw
       y1 = 0
       y2 = segh
   endif else begin
       x1 = detcols - ((i-9) * segw)
       x2 = x1 - segw
       y1 = detrows
       y2 = segh
   endelse
   detsec=string('[',x1,':',x2,',',y1,':',y2,']',format='(A,I0,A,I0,A,I0,A,I0,A)')
   SXADDPAR, seghdr, 'DETSEC', detsec, FORMAT='(A)',AFTER='DATASEC'
   SXADDPAR, seghdr, 'CCDSUM', '1 1', FORMAT='(A)', AFTER='DETSEC'    
   
   ; do LTV1 correctly here
   if (i LT 9) then begin
       ltv1 = string(0-(segw*(i-1)),format='(F10.3)')
   endif else begin
       ltv1 = string(detcols-((i-9)*segw),format='(F10.3)')
   endelse
   
   SXADDPAR, seghdr, 'LTV1', ltv1, FORMAT='(A)', AFTER='CCDSUM'
   
   if (i LT 9) then begin
       SXADDPAR, seghdr, 'LTV2', '-0.', FORMAT='(A)', AFTER='LTV1'
       SXADDPAR, seghdr, 'LTM1_1','1.', FORMAT='(A)', AFTER='LTV2'
       SXADDPAR, seghdr, 'LTM2_2','1.', FORMAT='(A)', AFTER='LTM1_1'
   endif else begin
       SXADDPAR, seghdr, 'LTV2','4044.', FORMAT='(A)', AFTER='LTV1'
       SXADDPAR, seghdr, 'LTM1_1','-1.', FORMAT='(A)', AFTER='LTV2'
       SXADDPAR, seghdr, 'LTM2_2','-1,', FORMAT='(A)', AFTER='LTM1_1'
   endelse
   
   SXADDPAR, seghdr, 'DTV1', 0, FORMAT='(I)', AFTER='LTM2_2'
   SXADDPAR, seghdr, 'DTV2', 0, FORMAT='(I)', AFTER='DTV1'
   SXADDPAR, seghdr, 'DTM1_1', 1, FORMAT='(I)', AFTER='DTV2'
   SXADDPAR, seghdr, 'DTM2_2', 1, FORMAT='(I)', AFTER='DTM1_1'

   if keyword_set(CHAN_ID) then im[*] = uint(i)
   FITS_WRITE, fcb, im, seghdr
endfor

FITS_CLOSE, fcb

; if image array pixel value is specified, set the image array pixels to that
if (verbose) then print,'Setting image array pixel values'
IF N_Elements(imval) EQ 1 THEN for i = 1, 16 do $
   seg_imval,filename,i,uint(imval),/image,PRE_COLS=pre_cols, $
      OVER_COLS=over_cols, OVER_ROWS=over_rows

; if overscan pixel value is specified, set the overscan pixels to that
if (verbose) then print,'Setting overscan pixel values'
IF N_Elements(overval) EQ 1 THEN for i = 1, 16 do $
   seg_overval, filename,i, overval, $
      PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
      VERBOSE=verbose

end