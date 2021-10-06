module modSlaterIndex
  implicit none
  private
  public :: SlaterIndex, SlaterMaximumCoincidence

  !TODO: rewrite to match complete OOP standard
  type SlaterIndex
    integer, dimension(:,:), allocatable  ::  spinOrbitals
    integer, dimension(:,:), allocatable  ::  spaceOrbitals
    real,    dimension(:),   allocatable  ::  coefs
    integer                               ::  NConf
  end type SlaterIndex

  interface SlaterIndex
    procedure ::  initSlater
  end interface SlaterIndex
  
  contains
  
  type(SlaterIndex) function initSlater(spinOrbitals)
    integer, dimension(:,:), intent(in) ::  spinOrbitals
    integer                           ::  i, j, Ne, NConf
    Ne = size(spinOrbitals,1)
    NConf = size(spinOrbitals,2)
    initSlater%NConf = NConf
    allocate(initSlater%spinOrbitals(Ne,NConf))
    allocate(initSlater%spaceOrbitals(Ne,NConf))
    
    initSlater%spinOrbitals = spinOrbitals
    do i=1,Ne
      do j=1,NConf
        initSlater%spaceOrbitals(i,j) = ceiling(spinOrbitals(i,j)/2.0)
      end do
    end do

    allocate(initSlater%coefs(NConf))
    !> Szabo Table 2.7
    !> TODO: Make sure the doubles are fed in the correct order
    if ( NConf == 2 ) then
      initSlater%coefs = 1.0/sqrt(2.0)
    elseif ( NConf == 6 ) then
      initSlater%coefs(1:2) =  2.0/sqrt(12.0)
      initSlater%coefs(3)   = -1.0/sqrt(12.0)
      initSlater%coefs(4:5) =  1.0/sqrt(12.0)
      initSlater%coefs(6)   = -1.0/sqrt(12.0)
    elseif ( NConf == 4 ) then
      initSlater%coefs = 1.0/2.0
    endif
    


  end function initSlater


  !> Input two spinOrbital configurations
  !> output the sgn associated with putting them in maximum coincindence
  !> output the spin Orbitals which differ as spinOrbI1, spinOrbI2, spinOrbIp1, spinOrbIp2

  subroutine SlaterMaximumCoincidence(spinOrbI, spinOrbIp, sgn, &
                                      numDiffering, spinOrbI1, spinOrbI2, spinOrbIp1, spinOrbIp2)
    integer,           intent(in)               ::  spinOrbI(:), spinOrbIp(:)
    integer,           intent(out)              ::  sgn, numDiffering, spinOrbI1, spinOrbI2, spinOrbIp1, spinOrbIp2
    integer                                     ::  i, indexLocation, numMatching, Nele
    logical                                     ::  firstDiff
    integer, dimension(:), allocatable          ::  tmpSpin
  

    Nele = size(spinOrbI,1)

    allocate(tmpSpin(Nele))
    tmpSpin = 0

    ! scan through slaterIp to see by how many spin indices it differs from slaterI
    numMatching=0
    do i=1,Nele
      if ( any( spinOrbIp(i)==spinOrbI ) ) then ! We have one matching
        numMatching = numMatching + 1
      endif
    enddo
    ! We need to report which spin orbitals in slaterI and slaterIp don't match
    numDiffering = Nele - numMatching
   
    if(numMatching > (Nele-3)) then !! Only care if we differ by less than 2

      tmpSpin = spinOrbIp
      sgn = 1
      do i=1,Nele
        if ( any( spinOrbI(i) == tmpSpin ) ) then
          indexLocation = findloc(tmpSpin,spinOrbI(i),1)
          if(indexLocation .ne. i) then 
            ! Only need to swap the elements if the index locations are not equal
            tmpSpin(indexLocation) = tmpSpin(i)
            tmpSpin(i) = spinOrbI(i)
          
            sgn = sgn * (-1)
          endif
        endif
      enddo

      ! |I>  = | ...mn...>
      ! |Ip> = | ...pq...>

      !> Assign the spin Orbitals which are different
      tmpSpin = tmpSpin - spinOrbI 
      firstDiff = .true.
      do i=1, Nele
        if ((tmpSpin(i) .ne. 0)  .and. (firstDiff .eqv. .true.)) then

          spinOrbI1 = spinOrbI(i)
          ! find the spin orbital in slater Ip
          indexLocation = findloc(spinOrbIp, tmpSpin(i) + spinOrbI(i),1)
          spinOrbIp1 = spinOrbIp(indexLocation)

          firstDiff=.false.
        elseif ((tmpSpin(i) .ne. 0)  .and. (firstDiff .eqv. .false.)) then
          
          spinOrbI2 = spinOrbI(i)
          ! find the spin orbital in slater Ip
          indexLocation = findloc(spinOrbIp, tmpSpin(i) + spinOrbI(i),1)
          spinOrbIp2 = spinOrbIp(indexLocation)
           
        endif
      enddo

      print *, ">> has <<"
      print *, "spinOrbI ", spinOrbI
      print *, "spinOrbIp", spinOrbIp
      !print *, "spaceOrbI1: ", spaceOrbI1, "spaceOrbI2: ", spaceOrbI2
      !print *, "spaceOrbIp1:", spaceOrbIp1, "spaceOrbIp2:", spaceOrbIp2
      print *, ">> with sign = ", sgn, ", numDiffering =",numDiffering,"<<"
      
    else
      sgn=0
      spinOrbI1 = 0
      spinOrbI2 = 0
      spinOrbIp1 = 0
      spinOrbIp2 = 0

    endif

  end subroutine SlaterMaximumCoincidence


  

end module modSlaterIndex
