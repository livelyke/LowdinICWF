module params
  use modSlaterIndex

  implicit none
  private
  integer, public                         ::  m1=1836,m2=1836, nDim, Nele, NeleSpinOrb, &
                                              NSlater, &
                                              NeleSpatialOrb, &
                                              NnucOrb, kickDir
  real, public                            ::  dxn = 0.1, nBoxL = 0.5, nBoxR = 6.00, &
                                              mu_n, dt=0.1, dtImag=0.01, finalTime=1000, kappa
  real, allocatable, dimension(:), public ::  nAxis
  character(len=26),public                ::  CITruncation
  type(SlaterIndex), dimension(:), allocatable,public  ::  slaterIndices
  logical,                      public    ::  debug, electronicOnly, kick

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
    read(1,*) electronicOnly, RVal
    read(1,*) NnucOrb
    read(1,*) dt, finalTime
    read(1,*) kick
    read(1,*) kappa, kickDir

    if (electronicOnly .eqv. .true.) then
      NnucOrb = 1
    endif

    if (debug .eqv. .true.) then
      print *, "CITruncation    : ", CITruncation
      print *, "Nele            : ", Nele
      print *, "NeleSpatialOrb  : ", NeleSpatialOrb
      print *, "debug           : ", debug
      print *, "electronicOnly  : ", electronicOnly
      print *, "RVal            : ", RVal
      print *, "NnucOrb         : ", NnucOrb
      write(*,'(A,F5.3,A,F6.0)') "dt, finalTime    : ", dt, ", ", finalTime
      print *, "kick            : ", kick
      print *, "kappa, kickDir  : ", kappa, kickDir
    endif

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

      !do i=1,NSlater
      !  do j=1,slaterIndices(i)%NConf 
      !    do s1=(i+1),NSlater
      !      do k=1,slaterIndices(s1)%NConf
      !  
      !        write(*,'(A,I3,A,I1,A,I3,A,I1)') "I=", i, ".",j,"  Ip=",s1,".",k
      !        call SlaterMaximumCoincidence( &
      !              slaterIndices(i)%spinOrbitals(:,j), &
      !              slaterIndices(s1)%spinOrbitals(:,k), sgn, numDiffering, &
      !              spaceOrbI1, spaceOrbI2, spaceOrbIp1, spaceOrbIp2)
      !        print *, "------------"
      !      enddo
      !    enddo
      !  enddo
      !enddo
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
    elseif ( CITruncation=='doubles-singlet' ) then
      i=0
      !> Same to Same
      do s1 = 1, Nele, 2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          i = i + 1
        enddo
      enddo

      !> Same to Different
      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        do s2=s1+1,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            do s2Promotion=s1Promotion+1,NeleSpatialOrb
              i = i+1
            enddo
          enddo
        enddo
      enddo

      !> Different to Same
      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        !> Besides The neighboring one which comes from the same space orb.
        do s2=s1+3,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
              i = i+1
          enddo
        enddo
      enddo

      !> Different to Different A and b
      do s1=1,Nele/2
        !> Promote every other opposite spin electron in the determinant
        do s2=s1+1,Nele/2
          !> Write s Promotions in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            do s2Promotion = s1Promotion+1,NeleSpatialOrb
              i = i+2
            enddo
          enddo
        enddo
      enddo

      NSlater = 1 + i !> Accounts for both spin channels in single and double 
    elseif ( CITruncation=='singles-doubles-singlet' ) then
      i=0
      !> Singlet and Same to Same
      do s1 = 1, Nele, 2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          i = i + 2 !< singlet and same to same
        enddo
      enddo

      !> Same to Different
      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        do s2=s1+1,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            do s2Promotion=s1Promotion+1,NeleSpatialOrb
              i = i+1
            enddo
          enddo
        enddo
      enddo

      !> Different to Same
      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        !> Besides The neighboring one which comes from the same space orb.
        do s2=s1+3,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
              i = i+1
          enddo
        enddo
      enddo

      !> Different to Different A and b
      do s1=1,Nele/2
        !> Promote every other opposite spin electron in the determinant
        do s2=s1+1,Nele/2
          !> Write s Promotions in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            do s2Promotion = s1Promotion+1,NeleSpatialOrb
              i = i+2
            enddo
          enddo
        enddo
      enddo

      NSlater = 1 + i !> Accounts for both spin channels in single and double 
    endif

    if (debug) then
      write(*,'(A,I5)') "NSlaters =", NSlater
    endif

  end function countSlaters

  subroutine populateSlaters(slaterIndices, CITruncation)
    type(SlaterIndex), intent(inout)        ::  slaterIndices(:)
    character(len=26), intent(in)     ::  CITruncation
    integer   :: i, j, a, b, r, s, aSpin, bSpin, rSpin, sSpin
    integer,dimension(:,:), allocatable       ::  unexcitedSpin, unexcitedSpace, spinSet
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
    elseif ( CITruncation=='doubles-singlet' ) then
      allocate(spinSet(Nele, 1))
      spinSet = 0
      !> Do the double excitation from same space orbital to same space orbital
      do s1=1,Nele,2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          spinSet(:,1) = unexcitedSpin(:,1)
          
          spinSet(s1,1) = s1Promotion
          spinSet(s1+1,1) = s1Promotion+1
          
          i = i+1
          slaterIndices(i) = SlaterIndex(spinSet)
        enddo
      enddo
      deallocate(spinSet)
      
      !> Do the double excitations from same to different
      allocate(spinSet(Nele, 2))
      spinSet = 0

      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        do s2=s1+1,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            do s2Promotion=s1Promotion+1,NeleSpatialOrb
              spinSet(:,1) = unexcitedSpin(:,1)
              spinSet(:,2) = unexcitedSpin(:,1)
              
              !> Map from Space orbitals to spin
              spinSet(s1,1) = 2*s1Promotion - 1
              spinSet(s2,1) = 2*s2Promotion
              spinSet(s1,2) = 2*s2Promotion - 1
              spinSet(s2,2) = 2*s1Promotion
              i = i+1
              slaterIndices(i) = SlaterIndex(spinSet)
            enddo
          enddo
        enddo
      enddo

      !> From Different to same
      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        !> Besides The neighboring one which comes from the same space orb.
        do s2=s1+3,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            s2Promotion=s1Promotion
            spinSet(:,1) = unexcitedSpin(:,1)
            spinSet(:,2) = unexcitedSpin(:,1)
            
            !> Map from Space orbitals to spin
            spinSet(s1,1) = 2*s1Promotion - 1
            spinSet(s2,1) = 2*s2Promotion
            spinSet(s1+1,2) = 2*s2Promotion - 1
            spinSet(s2-1,2) = 2*s1Promotion
            i = i+1
            slaterIndices(i) = SlaterIndex(spinSet)
          enddo
        enddo
      enddo

      deallocate(spinSet)

      !> From different to different A      
      allocate(spinSet(Nele, 6))
      spinSet = 0
      !> do loops over a,b,r,s
      do a=1,Nele/2
        !> Promote every other opposite spin electron in the determinant
        do b=a+1,Nele/2
          !> Write s Promotions in terms of space Orbitals
          do r = (Nele/2+1),NeleSpatialOrb
            do s = r+1,NeleSpatialOrb
              do j=1,6
                spinSet(:,j) = unexcitedSpin(:,1)
              enddo
              !> ab->rs
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b-1
              rSpin = 2*r-1
              sSpin = 2*s-1
              
              spinSet(aSpin,1) = rSpin 
              spinSet(bSpin,1) = sSpin 
              
              !> \bar{ab}->\bar{rs}
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b
              rSpin = 2*r
              sSpin = 2*s
              
              spinSet(aSpin,2) = rSpin 
              spinSet(bSpin,2) = sSpin 
              
              !> \bar{a}b->\bar{s}r
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,3) = sSpin 
              spinSet(bSpin,3) = rSpin 
              
              !> \bar{a}b->\bar{r}s
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,4) = rSpin 
              spinSet(bSpin,4) = sSpin 
              
              !> a\bar{b}->r\bar{s}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,5) = rSpin 
              spinSet(bSpin,5) = sSpin 
              
              !> a\bar{b}->s\bar{r}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,6) = sSpin 
              spinSet(bSpin,6) = rSpin 
            
              i = i+1
              slaterIndices(i) = SlaterIndex(spinSet)
              
            enddo
          enddo
        enddo
      enddo
      
      deallocate(spinSet)

      !> From different to different B
      allocate(spinSet(Nele, 4))
      spinSet = 0
      !> do loops over a,b,r,s
      do a=1,Nele/2
        !> Promote every other opposite spin electron in the determinant
        do b=a+1,Nele/2
          !> Write s Promotions in terms of space Orbitals
          do r = (Nele/2+1),NeleSpatialOrb
            do s = r+1,NeleSpatialOrb
              do j=1,4
                spinSet(:,j) = unexcitedSpin(:,1)
              enddo
              
              !> \bar{a}b->\bar{s}r
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,1) = sSpin 
              spinSet(bSpin,1) = rSpin 
              
              !> \bar{a}b->\bar{r}s
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,2) = rSpin 
              spinSet(bSpin,2) = sSpin 
              
              !> a\bar{b}->r\bar{s}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,3) = rSpin 
              spinSet(bSpin,3) = sSpin 
              
              !> a\bar{b}->s\bar{r}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,4) = sSpin 
              spinSet(bSpin,4) = rSpin 
            
              i = i+1
              slaterIndices(i) = SlaterIndex(spinSet)
              
            enddo
          enddo
        enddo
      enddo
      
    elseif ( CITruncation=='singles-doubles-singlet' ) then
      !> Singles
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
      deallocate(spinSet)
 
      allocate(spinSet(Nele, 1))
      spinSet = 0
      !> Do the double excitation from same space orbital to same space orbital
      do s1=1,Nele,2
        do s1Promotion = (Nele+1),NeleSpinOrb,2
          spinSet(:,1) = unexcitedSpin(:,1)
          
          spinSet(s1,1) = s1Promotion
          spinSet(s1+1,1) = s1Promotion+1
          
          i = i+1
          slaterIndices(i) = SlaterIndex(spinSet)
        enddo
      enddo
      deallocate(spinSet)
      
      !> Do the double excitations from same to different
      allocate(spinSet(Nele, 2))
      spinSet = 0

      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        do s2=s1+1,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            do s2Promotion=s1Promotion+1,NeleSpatialOrb
              spinSet(:,1) = unexcitedSpin(:,1)
              spinSet(:,2) = unexcitedSpin(:,1)
              
              !> Map from Space orbitals to spin
              spinSet(s1,1) = 2*s1Promotion - 1
              spinSet(s2,1) = 2*s2Promotion
              spinSet(s1,2) = 2*s2Promotion - 1
              spinSet(s2,2) = 2*s1Promotion
              i = i+1
              slaterIndices(i) = SlaterIndex(spinSet)
            enddo
          enddo
        enddo
      enddo

      !> From Different to same
      do s1=1,Nele,2
        !> Promote every other opposite spin electron in the determinant
        !> Besides The neighboring one which comes from the same space orb.
        do s2=s1+3,Nele,2
          !> Write s1Promotion in terms of space Orbitals
          do s1Promotion = (Nele/2+1),NeleSpatialOrb
            s2Promotion=s1Promotion
            spinSet(:,1) = unexcitedSpin(:,1)
            spinSet(:,2) = unexcitedSpin(:,1)
            
            !> Map from Space orbitals to spin
            spinSet(s1,1) = 2*s1Promotion - 1
            spinSet(s2,1) = 2*s2Promotion
            spinSet(s1+1,2) = 2*s2Promotion - 1
            spinSet(s2-1,2) = 2*s1Promotion
            i = i+1
            slaterIndices(i) = SlaterIndex(spinSet)
          enddo
        enddo
      enddo

      deallocate(spinSet)

      !> From different to different A      
      allocate(spinSet(Nele, 6))
      spinSet = 0
      !> do loops over a,b,r,s
      do a=1,Nele/2
        !> Promote every other opposite spin electron in the determinant
        do b=a+1,Nele/2
          !> Write s Promotions in terms of space Orbitals
          do r = (Nele/2+1),NeleSpatialOrb
            do s = r+1,NeleSpatialOrb
              do j=1,6
                spinSet(:,j) = unexcitedSpin(:,1)
              enddo
              !> ab->rs
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b-1
              rSpin = 2*r-1
              sSpin = 2*s-1
              
              spinSet(aSpin,1) = rSpin 
              spinSet(bSpin,1) = sSpin 
              
              !> \bar{ab}->\bar{rs}
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b
              rSpin = 2*r
              sSpin = 2*s
              
              spinSet(aSpin,2) = rSpin 
              spinSet(bSpin,2) = sSpin 
              
              !> \bar{a}b->\bar{s}r
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,3) = sSpin 
              spinSet(bSpin,3) = rSpin 
              
              !> \bar{a}b->\bar{r}s
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,4) = rSpin 
              spinSet(bSpin,4) = sSpin 
              
              !> a\bar{b}->r\bar{s}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,5) = rSpin 
              spinSet(bSpin,5) = sSpin 
              
              !> a\bar{b}->s\bar{r}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,6) = sSpin 
              spinSet(bSpin,6) = rSpin 
            
              i = i+1
              slaterIndices(i) = SlaterIndex(spinSet)
              
            enddo
          enddo
        enddo
      enddo
      
      deallocate(spinSet)

      !> From different to different B
      allocate(spinSet(Nele, 4))
      spinSet = 0
      !> do loops over a,b,r,s
      do a=1,Nele/2
        !> Promote every other opposite spin electron in the determinant
        do b=a+1,Nele/2
          !> Write s Promotions in terms of space Orbitals
          do r = (Nele/2+1),NeleSpatialOrb
            do s = r+1,NeleSpatialOrb
              do j=1,4
                spinSet(:,j) = unexcitedSpin(:,1)
              enddo
              
              !> \bar{a}b->\bar{s}r
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,1) = sSpin 
              spinSet(bSpin,1) = rSpin 
              
              !> \bar{a}b->\bar{r}s
              !> map from space to spin
              aSpin = 2*a
              bSpin = 2*b-1
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,2) = rSpin 
              spinSet(bSpin,2) = sSpin 
              
              !> a\bar{b}->r\bar{s}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r-1
              sSpin = 2*s
              
              spinSet(aSpin,3) = rSpin 
              spinSet(bSpin,3) = sSpin 
              
              !> a\bar{b}->s\bar{r}
              !> map from space to spin
              aSpin = 2*a-1
              bSpin = 2*b
              rSpin = 2*r
              sSpin = 2*s-1
              
              spinSet(aSpin,4) = sSpin 
              spinSet(bSpin,4) = rSpin 
            
              i = i+1
              slaterIndices(i) = SlaterIndex(spinSet)
              
            enddo
          enddo
        enddo
      enddo
    
    endif 

  end subroutine populateSlaters


end module params
