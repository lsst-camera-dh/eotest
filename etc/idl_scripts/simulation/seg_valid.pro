
function seg_valid, infile, exten, WRITE=write, MESSAGE=message

compile_opt IDL2

  ; check if the input file exists and is readable, if not exit
  if (~file_test(infile, /READ)) then begin
    if keyword_set(MESSAGE) then message='Error: Specified file does not exist or is unreadable'
    return, -1
  endif
  
  ; check if the input file exists and is writable, if not exit
  if keyword_set(WRITE) then begin
     if (~file_test(infile, /WRITE)) then begin
        if keyword_set(MESSAGE) then message='Error: Specified file does not exist or is unwritable'
        return, -1
     endif
  endif
  
  ; check if the segment requested makes sense
  if ((exten LT 1) || (exten GT 16)) then begin
    if keyword_set(MESSAGE) then message='Error: Specified extension number is invalid'
    return, -1
  endif
  
  ; if everything checked out, return 0
  return, 0

end