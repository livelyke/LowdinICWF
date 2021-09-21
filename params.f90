module params
  use modSlaterIndex

  implicit none
  private
  integer, public                         ::  m1=1836,m2=1836, nDim, Nele=4, NeleSpinOrb, &
                                              NSlater, &
                                              NeleSpatialOrb=4, &
                                              NnucOrb=1
  real, public                            ::  dxn = 0.1, nBoxL = 0.5, nBoxR = 10, mu_n, dt=0.01, dtImag=0.1
  real, allocatable, dimension(:), public ::  nAxis
  character(len=20),public                ::  CITruncation='doubles'
  type(SlaterIndex), dimension(:), allocatable,public  ::  slaterIndices

  public :: assignConsts
  contains
  
  subroutine assignConsts()
    integer   :: i
    integer,dimension(:), allocatable       ::  unexcited, set
    integer                                 ::  s1, s1Promotion, s2, s2Promotion
    
    !integer               :: spaceOrbI1, spaceOrbI2, spaceOrbIp1, spaceOrbIp2,sgn, numDiffering
    mu_n=m1*m2/(m1+m2)
    
    nDim = (nBoxR-nBoxL)/dxn
    allocate(nAxis(nDim))
    nAxis = (/ (nBoxL + i*dxn, i=0,nDim-1) /)

    NeleSpinOrb = NeleSpatialOrb*2

    !TODO: Add Spin adapted configurations

    if ( CITruncation == 'doubles' ) then
    !TODO:  This is only accurate for moving an up and a down spin
    !       Needs to be expanded for double up/down promotion above Nele
    !       double spin down/up promotion
    !       I'm to lazy to do the combinatorics, just count how many 
    !       unique ones there are for a given spin

      ! This counting procedure excludes `identical' bases
      ! e.g. for Nele=4, NeleSpinOrb=8
      ! included would be 5274
      ! excluded would be 7254
      ! because 7254 = -5274
      ! and this doesn't contribute to the span of the basis.
      i=0
      do s1 = 1, Nele, 2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          do s2 = s1+2, Nele, 2
            do s2Promotion = (s1Promotion+2), NeleSpinOrb,2
              i = i+1
            enddo
          enddo
        enddo
      enddo 
    
      NSlater = 1 + (Nele/2 * (NeleSpinOrb-Nele)/2)**2 &
               + 2*i ! 
      
      allocate(slaterIndices(NSlater))
      allocate(unexcited(Nele))
      unexcited = (/ (i, i=1,Nele) /)
      allocate(set(Nele))
 
      ! do first all the double excitations involving the promotion
      ! of an up spin and a down spin
      i=1
      slaterIndices(i) = SlaterIndex(unexcited)
      do s1 = 1, Nele, 2 
        ! start from lowest spin orbital
        ! promote to every possible spin orbital of same spin
        do s1Promotion = (Nele+1), NeleSpinOrb, 2
          do s2 = 2, Nele, 2
            do s2Promotion = (Nele+2), NeleSpinOrb, 2
              set = unexcited
              set(s1) = s1Promotion
              set(s2) = s2Promotion
              i = i+1
              slaterIndices(i) = SlaterIndex(set)
            end do
          end do
        end do
      end do

      ! do now all the excitations involving the promotion 
      ! of two up or two down spins
      ! first for 'up' spins
      do s1 = 1, Nele, 2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          do s2 = s1+2, Nele, 2
            do s2Promotion = (s1Promotion+2), NeleSpinOrb,2
              set = unexcited
              set(s1) = s1Promotion
              set(s2) = s2Promotion
              i = i+1
              slaterIndices(i) = SlaterIndex(set)
            enddo
          enddo
        enddo
      enddo 
      
      ! now for 'down' spins
      do s1 = 2, Nele, 2
        do s1Promotion = (Nele+2),NeleSpinOrb,2
          do s2 = s1+2, Nele, 2
            do s2Promotion = (s1Promotion+2), NeleSpinOrb,2
              set = unexcited
              set(s1) = s1Promotion
              set(s2) = s2Promotion
              i = i+1
              slaterIndices(i) = SlaterIndex(set)
            enddo
          enddo
        enddo
      enddo 
    end if

    ! Sep 17, debug the reporting of the differing spinOrbitals
    !do i=1,NSlater
    !  do s1=(i+1),NSlater
    !    call SlaterMaximumCoincidence(slaterIndices,i, s1, Nele, sgn, numDiffering, &
    !                                  spaceOrbI1, spaceOrbI2, spaceOrbIp1, spaceOrbIp2)
    !  enddo
    !enddo


  end subroutine assignConsts
 !nAxis = (nBoxL:dxn:nBoxR); nDim = max(size(nAxis));

!m1 = 1836; m2 = 1836; mu_n = m1*m2/(m1+m2);


end module params
