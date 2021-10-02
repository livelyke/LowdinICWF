module params
  use modSlaterIndex

  implicit none
  private
  integer, public                         ::  m1=1836,m2=1836, nDim, Nele, NeleSpinOrb, &
                                              NSlater, &
                                              NeleSpatialOrb, &
                                              NnucOrb=1
  real, public                            ::  dxn = 0.1, nBoxL = 0.5, nBoxR = 10, mu_n, dt=0.01, dtImag=0.01
  real, allocatable, dimension(:), public ::  nAxis
  character(len=26),public                ::  CITruncation
  type(SlaterIndex), dimension(:), allocatable,public  ::  slaterIndices
  logical,                      public    ::  debug, electronicOnly

  public :: assignConsts
  contains

  subroutine readInp(CITruncation,Nele,NeleSpatialOrb,debug,electronicOnly,RVal)
    character(len=26), intent(inout)  ::  CITruncation
    integer,           intent(inout)  ::  Nele, NeleSpatialOrb
    logical,           intent(inout)  ::  debug, electronicOnly
    real,              intent(out)    ::  RVal

    open(1,file='inpICWF',status='old')
    read(1,*) CITruncation
    read(1,*) Nele
    read(1,*) NeleSpatialOrb
    read(1,*) debug
    read(1,*) electronicOnly
    if(electronicOnly .eqv. .true.) then
      read(1,*) RVal
    endif

    print *, "CITruncation ", CITruncation
    print *, "Nele ", Nele
    print *, "NeleSpatialOrb ", NeleSpatialOrb
    print *, "debug ", debug
    print *, "electronicOnly ", electronicOnly
    print *, "RVal ", RVal

  end subroutine readInp
  
  subroutine assignConsts()
    integer               :: i 
    integer               :: spaceOrbI1, spaceOrbI2, spaceOrbIp1, spaceOrbIp2, &
                             s1, sgn, numDiffering
    real                  :: Rval
    mu_n=m1*m2/(m1+m2)
   

    call readInp(CITruncation,Nele,NeleSpatialOrb,debug,electronicOnly,RVal)


    if (.not.electronicOnly) then 
      nDim = (nBoxR-nBoxL)/dxn
      allocate(nAxis(nDim))
      nAxis = (/ (nBoxL + i*dxn, i=0,nDim-1) /)
    else
      nDim=1
      allocate(nAxis(nDim))
      nAxis(1) = RVal
      close(1)
      NnucOrb = 1
    endif

    NeleSpinOrb = 2*NeleSpatialOrb

    !TODO: Add Spin adapted configurations
    !TODO: Add singles

    NSlater = countSlaters(CITruncation)
    allocate(slaterIndices(NSlater))

    call populateSlaters(slaterIndices, CITruncation)

    if (debug .eqv. .true.) then
      do i=1,NSlater
        print *, "------------"
        print *, "Slater I=", i
        print *, "spinOrbs ", slaterIndices(i)%spinOrbitals
        print *, "spaceOrbs", slaterIndices(i)%spaceOrbitals
      enddo
      print *, "------------"

      do i=1,NSlater
        do s1=(i+1),NSlater
          call SlaterMaximumCoincidence(slaterIndices,i, s1, Nele, sgn, numDiffering, &
                                        spaceOrbI1, spaceOrbI2, spaceOrbIp1, spaceOrbIp2)
        enddo
      enddo
    endif

  end subroutine assignConsts

  function countSlaters(CITruncation) result(NSlater)

    character(len=26), intent(in)     ::  CITruncation
    integer                   ::  NSlater
    integer                   ::  i, s1, s1Promotion, s2, s2Promotion
    

    if ( CITruncation .eq. 'doubles-singlet' ) then
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
    endif
    if ( CITruncation .eq. 'singles-doubles-singlet' ) then
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
      do s1 = 1, Nele, 2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          i = i + 1
        enddo
      enddo 
      
      NSlater = 1 + (Nele/2 * (NeleSpinOrb-Nele)/2)**2 &
               + 2*i !> Accounts for both spin channels in single and double 
    endif

  end function countSlaters

  subroutine populateSlaters(slaterIndices,CITruncation)
    type(SlaterIndex), intent(inout)        ::  slaterIndices(:)
    character(len=26), intent(in)     ::  CITruncation
    integer   :: i
    integer,dimension(:), allocatable       ::  unexcited, set
    integer                                 ::  s1, s1Promotion, s2, s2Promotion
 
    !allocate(slaterIndices(NSlater))
    allocate(unexcited(Nele))
    unexcited = (/ (i, i=1,Nele) /)
    allocate(set(Nele))
 
    ! do first all the double excitations involving the promotion
    ! of an up spin and a down spin
    i=1
    slaterIndices(i) = SlaterIndex(unexcited)

    if (CITruncation=='doubles-singlet') then
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
    elseif (CITruncation=='singles-doubles-singlet') then
      !> Start with singles excitation
      do s1 = 1, Nele, 2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          set = unexcited
          set(s1) = s1Promotion
          i = i+1
          slaterIndices(i) = SlaterIndex(set)
        enddo
      enddo 
      !> Now for opposite spin channel
      do s1 = 2, Nele, 2
        do s1Promotion = (Nele+2),NeleSpinOrb,2
          set = unexcited
          set(s1) = s1Promotion
          i = i+1
          slaterIndices(i) = SlaterIndex(set)
        enddo
      enddo 

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

    endif 

  end subroutine populateSlaters


end module params
