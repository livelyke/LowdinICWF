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
    integer               :: i, j, k
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
        print *, "spinOrbs " 
        do s1=1,size(slaterIndices(i)%spinOrbitals,2)
          print *, slaterIndices(i)%spinOrbitals(:,s1)
        enddo
        print *, "spaceOrbs"
        do s1=1,size(slaterIndices(i)%spaceOrbitals,2)
          print *, slaterIndices(i)%spaceOrbitals(:,s1)
        enddo
      enddo
      print *, "------------"

      do i=1,NSlater
        do j=1,slaterIndices(i)%NConf 
          do s1=(i+1),NSlater
            do k=1,slaterIndices(s1)%NConf

              print *, "SlaterI  =", i, "conf =",j 
              print *, "SlaterIp =", s1, "conf =",k 
              call SlaterMaximumCoincidence( &
                    slaterIndices(i)%spinOrbitals(:,j), &
                    slaterIndices(s1)%spinOrbitals(:,k), sgn, numDiffering, &
                    spaceOrbI1, spaceOrbI2, spaceOrbIp1, spaceOrbIp2)
              print *, "------------"
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine assignConsts

  function countSlaters(CITruncation) result(NSlater)

    character(len=26), intent(in)     ::  CITruncation
    integer                   ::  NSlater
    integer                   ::  i, s1, s1Promotion, s2, s2Promotion
    
    if ( CITruncation .eq. 'singles-singlet' ) then
      i=0
      do s1 = 1, Nele, 2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          i = i + 1
        enddo
      enddo 
      NSlater = 1 + i !> Accounts for both spin channels in single and double 
    endif

  end function countSlaters

  subroutine populateSlaters(slaterIndices, CITruncation)
    type(SlaterIndex), intent(inout)        ::  slaterIndices(:)
    character(len=26), intent(in)     ::  CITruncation
    integer   :: i
    integer,dimension(:,:), allocatable       ::  unexcitedSpin, unexcitedSpace, spaceSet, spinSet
    integer                                 ::  s1, s1Promotion, s2, s2Promotion
 
    !allocate(slaterIndices(NSlater))
    allocate(unexcitedSpin(Nele,1))
    unexcitedSpin(:,1) = (/ (i, i=1,Nele) /)
    allocate(unexcitedSpace(Nele,1))
    unexcitedSpace(:,1) = (/ (ceiling(i/2.0), i=1,Nele) /)
 
    ! do first all the double excitations involving the promotion
    ! of an up spin and a down spin
    i=1
    slaterIndices(i) = SlaterIndex(unexcitedSpin)

    if ( CITruncation=='singles-singlet') then
      allocate(spaceSet(Nele, 2))
      spaceSet = 0
      allocate(spinSet(Nele, 2))
      spinSet = 0
      !> go through spin orbitals in pairs
      do s1 = 1, Nele ,2
        do s1Promotion = (Nele+1), NeleSpinOrb, 2
          spinSet(:,1) = unexcitedSpin(:,1)
          spinSet(:,2) = unexcitedSpin(:,1)

          !> spin up and spin down channels
          do s2=0,1
            spinSet(s1+s2,1+s2) = s1Promotion+s2
          enddo 
          i = i+1
          slaterIndices(i) = SlaterIndex(spinSet)
        enddo
      enddo 
 
    endif 

  end subroutine populateSlaters


end module params
