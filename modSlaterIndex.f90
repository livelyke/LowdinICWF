module modSlaterIndex
  implicit none
  private
  public :: SlaterIndex, SlaterMaximumCoincidence

  !TODO: rewrite to match complete OOP standard
  type SlaterIndex
    integer, dimension(:), allocatable  ::  spinOrbitals
    integer, dimension(:), allocatable  ::  spaceOrbitals
  end type SlaterIndex

  interface SlaterIndex
    procedure ::  initSlater
  end interface SlaterIndex
  
  contains
  
  type(SlaterIndex) function initSlater(spinOrbitals)
    integer, dimension(:), intent(in) ::  spinOrbitals
    integer                           ::  i, Ne
    Ne = size(spinOrbitals)
    allocate(initSlater%spinOrbitals(Ne))
    allocate(initSlater%spaceOrbitals(Ne))
    
    initSlater%spinOrbitals = spinOrbitals
    do i=1,Ne
      initSlater%spaceOrbitals(i) = ceiling(spinOrbitals(i)/2.0)
    end do

  end function initSlater

  subroutine SlaterMaximumCoincidence(slaterIndices, slaterI, slaterIp, Nele, sgn, &
                                      numDiffering, spinOrbI1, spinOrbI2, spinOrbIp1, spinOrbIp2)
    type(SlaterIndex),   intent(in)             ::  slaterIndices(:)
    integer,           intent(in)               ::  Nele, slaterI, slaterIp
    integer,           intent(out)              ::  sgn, numDiffering, spinOrbI1, spinOrbI2, spinOrbIp1, spinOrbIp2
    integer                                     ::  i, indexLocation, numMatching
    logical                                     ::  firstDiff
    integer, dimension(:), allocatable          ::  spaceOrbIp, spinOrbIp, spaceOrbI, spinOrbI, tmpSpin

    allocate(spinOrbI(Nele))
    !allocate(spaceOrbI(Nele))
    allocate(spinOrbIp(Nele))
    !allocate(spaceOrbIp(Nele))
    allocate(tmpSpin(Nele))
    spinOrbI    =   slaterIndices(slaterI)%spinOrbitals;
    !spaceOrbI   =   slaterIndices(slaterI)%spaceOrbitals;
    spinOrbIp   =   slaterIndices(slaterIp)%spinOrbitals;
    !spaceOrbIp  =   slaterIndices(slaterIp)%spaceOrbitals;


    ! scan through slaterIp to see by how many spin indices it differs from slaterI
    numMatching=0
    do i=1,Nele
      if ( any( spinOrbIp(i)==spinOrbI ) ) then ! We have one matching
        numMatching = numMatching + 1
      endif
    enddo
   
    ! The structure of the promotions may help with how to permute indices 
    ! Recall for spin up and down promotions
    !do s1 = 1, Nele, 2 
      !do s1Promotion = (Nele+1), NeleSpinOrb, 2
        !do s2 = 2, Nele, 2
          !do s2Promotion = (Nele+2), NeleSpinOrb, 2

    ! can also just sort arrays from smallest to largest . . .
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

      ! Okay, now they've been put into maximum coincidence, we need to
      ! report which spin orbitals in slaterI and slaterIp don't match
      numDiffering = Nele - numMatching

      ! |I>  = | ...mn...>
      ! |Ip> = | ...pq...>

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

      !print *, "-----------"
      !print *, "spinOrbI ", spinOrbI
      !print *, "spinOrbIp", tmpSpin + spinOrbI
      !print *, "spaceOrbI1: ", spaceOrbI1, "spaceOrbI2: ", spaceOrbI2
      !print *, "spaceOrbIp1:", spaceOrbIp1, "spaceOrbIp2:", spaceOrbIp2
      
    else
      sgn=0
      spinOrbI1 = 0
      spinOrbI2 = 0
      spinOrbIp1 = 0
      spinOrbIp2 = 0

    endif

  end subroutine SlaterMaximumCoincidence

  

end module modSlaterIndex
