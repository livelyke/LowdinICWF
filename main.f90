program icwf
  use params, only : assignConsts
  use operators, only : setOperators, H
  use tdEvolve, only : C, tdEvolveImag

  implicit none

  call assignConsts()

  call setOperators()

  !call tdEvolveImag(H) 
 
end program icwf
