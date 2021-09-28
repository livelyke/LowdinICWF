module tdEvolve
  use params,     only : NSlater,NnucOrb
  use operators,  only : H
  implicit none
  private
  public  ::  tdEvolveImag

  real, public, allocatable, dimension(:) :: C

  contains 

  subroutine tdEvolveImag(H)
    real, intent(in)  ::  H(:,:)
    integer           ::  i

    allocate(C(NSlater*NnucOrb))
    call random_number(C)
    !C = (/ (rand(2*i), i=1,NSlater*NnucOrb) /)
    C = C/sqrt(sum(C*C))
    print *, C
    print *, sqrt(sum(C*C))

  end subroutine 

end module tdEvolve
