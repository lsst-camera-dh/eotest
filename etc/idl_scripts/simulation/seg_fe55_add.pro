
pro seg_fe55_add, infile, segment, $
   EVENTS=events, BETAS=betas, GAIN=gain, $
   PRE_COLS=pre_cols, OVER_COLS=over_cols, OVER_ROWS=over_rows, $
   VERBOSE=verbose

;+
; NAME:
;   SEG_FE55_ADD
;
; PURPOSE:
;   Add simulated fe55 x-ray charge packets to an image segment
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;    seg_fe55_add, infile, segment, value [ , VERBOSE=verbose ]
;
; INPUTS:
;    INFILE:  filename of the multi extension FITS file to operate on
;    SEGMENT: number of the extension to operate on
;
; OPTIONAL INPUTS:

;
; INPUT KEYWORD PARAMETERS:
;    EVENTS : the total number of x-ray charge packets to add, default=10000
;    BETAS  : fraction of events that should be betas, nature=12% = 0.12
;    GAIN   ; system gain in e-/DN, defualt=1.00
;    PRE_COLS  : number of prescan columns in image segment, default=10
;    OVER_COLS : number of overscan columns in image segment, default=20
;    OVER_ROWS : number of overscan rows in image segment, default=20
;    VERBOSE   : if set, print some diagnostics
;
; OUTPUT KEYWORD PARAMETERS:
;
; EXAMPLES:
;
; WARNING:
;
; PROCEDURES USED:
;
; HISTORY:
; Written by: P. Doherty  April, 2013
;
;-
;-----------------------------------------------------------------------------
compile_opt IDL2

if N_Params() NE 2 then begin
  print,'SEG_FE55_ADD: Incorrect number of parameters. Exiting.'
  return
endif

if N_Elements(events)    EQ 0 then events    = 10000   ; number of x-ray events
if N_Elements(betas)     EQ 0 then betas     = 0       ; fraction of beta events to include (~12% in real world)
if N_Elements(gain)      EQ 0 then gain      = 1.00    ; system gain e-/DN
IF N_Elements(pre_cols)  EQ 0 THEN pre_cols  = 10
IF N_Elements(over_cols) EQ 0 THEN over_cols = 20
IF N_Elements(over_rows) EQ 0 THEN over_rows = 20

errmsg='No Error.'
if (seg_valid(infile, segment, message=errmsg) NE 0) then begin
  print,errmsg
  return
end

alphas = events * (1 - betas)
betas  = events * (betas)
; print,alphas,betas

fits_read, infile, im, hdr, exten_no=segment
im = long(im) ; temporarily make image a long to do the math without clipping
s = size(im)

; create an image array to hold the calculated signal pixels
; (just the image aray pixels)
xim = im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] 
xim = float(xim)
xim[*] = 0.00
sx = size(xim)

if (alphas GE 1) then begin
   print, 'Generating ', alphas, ' alpha events'
   for i = 1, alphas do begin
      x = uint(randomu(seed) * (sx[1])) 
      y = uint(randomu(seed) * (sx[2]))
      xim[x,y] = xim[x,y] + (1620 / gain)
   endfor
endif

if (betas GE 1) then begin
   print, 'Generating ', betas, ' beta events'
   for i = 1, betas do begin
      x = uint(randomu(seed) * (sx[1]))
      y = uint(randomu(seed) * (sx[2]))
      xim[x,y] = xim[x,y] + (1778 / gain)
   endfor
endif

; add fe55 signal image into image array of input image
im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] = $
    im[pre_cols:(s[1]-(over_cols+1)), 0:(s[2]-(over_rows+1))] + xim

; clip the result
im = im < 65535
im = im > 0

; convert back to a UINT
im = uint(im)

; save the file
MODFITS, infile, im, hdr, EXTEN_NO=segment

end

