program icwf
  use params, only : assignConsts, finalTime
  use operators, only : setOperators, H
  use tdEvolve, only : tdEvolveImag, tdEvolveReal, C

  implicit none
  
  call assignConsts()

  call setOperators()

  call tdEvolveImag(H)
 
  if (finalTime .ne. 0) then 
    call tdEvolveReal(H,C) 
  endif
 
end program icwf
