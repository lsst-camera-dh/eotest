
function bright_pix, file, PCT=pct, EXT=ext
    
if ~keyword_set(EXT) then ext=1
if ~keyword_set(PCT) then pct=90.0

fits_read, file, im, hdr, EXTEN_NO=ext
im    = float(im[10:521, 0:2002])  ; extract image area
bmean = mean(float(im[522:*, *]))  ; calculate bias mean
image = image - bmean              ; subtract bias mean from image
col_means=fltarr(512)              ; storage for column mean values
    
; find mean value of each CCD column
for i = 0, 511 do col_means[i] = mean(image[i,*])
    
; find any columns with mean > PCT of image mean
bcols = where(col_means GE (mean(image) * pct), count)


    
    bright_pix = where(((*info.imptr)) GT threshold, bp_count)
    if (bp_count NE 0) then bp_mean = mean((*info.imptr)[bright_pix]) else bp_mean = 0.0

    ; find bright columns
    colavgs = fltarr(imsize[1])
    for i = 0, imsize[1]-1 do colavgs[i] = median((*info.imptr)[i,*])
    bright_cols = where(colavgs GT threshold, bc_count)

    ; weed out the bright pixels that are also in bright columns or rows
    if (bc_count GT 0) then begin
        for i = 0, bc_count-1 do begin
            cols = bright_pix mod imsize[1]
            index = where(cols NE bright_cols[i], count)
            if (count NE 0) then bright_pix = bright_pix[index] else bright_pix = -1L
        endfor
    endif

    ; recalulate the number of defects
    if ((n_elements(bright_pix) EQ 1) AND (bright_pix[0] EQ -1)) then $
        bp_count = 0 else bp_count = n_elements(bright_pix)


    ; recalulate the number of defects
    if ((n_elements(bright_pix) EQ 1) AND (bright_pix[0] EQ -1)) then $
        bp_count = 0 else bp_count = n_elements(bright_pix)

    ; compute the mean value of the defects
    if (bp_count NE 0) then bp_mean = mean((*info.imptr)[bright_pix]) else bp_mean = 0.0

    return,bright_pix

end