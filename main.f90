program icwf
  use params, only : assignConsts
  use operators, only : setOperators, H
  use tdEvolve, only : tdEvolveImag

  implicit none
  real, allocatable, dimension(:,:) :: KDx, KDy, KDz
  
  call assignConsts()

  call setOperators()

  !call tdEvolveImag(H) 
 
end program icwf
