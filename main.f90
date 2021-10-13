program icwf
  use params, only : assignConsts
  use operators, only : setOperators, H
  use tdEvolve, only : tdEvolveImag, tdEvolveReal, C

  implicit none
  
  call assignConsts()

  call setOperators()

  call tdEvolveImag(H)
  
  call tdEvolveReal(H,C) 
 
end program icwf
