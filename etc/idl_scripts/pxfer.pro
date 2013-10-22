pro pxfer, expmin, expmax, step, filebase, EXT=ext

if ~keyword_set(EXT) then ext=1

; results are an array: [bias_mean, bias_rms, flat_mean, flat_var, gain, noise]
result = fltarr(6)          ; the result of a single measurement
expcount = 1 + (expmax - expmin)/step
print,expcount
data = fltarr(n_elements(result),expcount)   ; the complete result data set

outfile = string(filebase,'_',ext,'.txt', format='(A,A,I02,A)')
openw, ofile, outfile, /get_lun

printf,ofile,' Exp     Bias Mean   Bias RMS    Flat Mean   Flat Var   Gain e-/DN   Noise e-  '
printf,ofile,' ------ ----------- ----------- ----------- ----------- ----------- ----------- '

exp = 0
i=0
while (exp LE expmax) do begin
    flat = string(filebase,'_',exp,'s', format='(A,A,F06.2,A)') 
    bias = string(filebase, format='(A)')            
    pair_stats, flat, ext, result
    print,string(exp,result[0],result[1],result[2],result[3],result[4],result[5],format='(F5.2,6(F12.3))')
    printf,ofile,string(exp,result[0],result[1],result[2],result[3],result[4],result[5],format='(F5.2,6(F12.3))')
    data[*,i] = result  
    exp=exp+step
    i=i+1
endwhile

print,'Finished writing output file : ',outfile
close, ofile
free_lun, ofile

end

