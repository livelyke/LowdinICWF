module operators
  use params, only          : mu_n, dxn, nDim, Nele, NeleSpinOrb, NSlater, NeleSpatialOrb, &
                              NnucOrb, nAxis, CITruncation, slaterIndices
  use bases, only           : diagSymMatrix
  use modSlaterIndex, only  : SlaterIndex

  implicit none (type, external)
  external  :: ssyev, sgemm, sgeev, ssyevr, amux
  private
  public    :: setOperators

  real, dimension(:,:), allocatable, public   ::  H

  contains

  subroutine setOperators()
    real, dimension(:), allocatable             :: eigenVals, work, iwork
    integer, dimension(:), allocatable          :: isuppz, spaceIndI, spinIndI, nucRange, &
                                                   slaterIJ, slaterIpJp, spaceIndIp, spinIndIp
    real, dimension(:,:), allocatable           ::  nLapl, Tn, chi, Vnn, KTn, KVnn, &
                                                    IdN, KTe, KVee, MTn, MTe, MVee, &
                                                    MVen, MVnn, tmpMat
    real, dimension(:,:,:), allocatable         ::  KVen
    real, dimension(:), allocatable     ::  tmpVec
    real                                ::  MTeVal=0, MVeeVal =0, rval
    integer                             ::  i, j, k, l, m, n, o, p, lwork, info, slaterIpIp, &
                                            NIJ, il, iu, vl, vu, Ri
    integer                             ::  ii, ij, ji, jj, slaterI, slaterIp, slaterJ, slaterJp
    integer                             ::  TeFile, VeeFile, VenFile, sizeFile
    logical                             ::  continueReading = .true.
    complex                             ::  zval
    character(len=256)           ::  filepath

    allocate(nLapl(nDim,nDim))
    nLapl = 0
    allocate(Tn(nDim,nDim))
    Tn = 0
    allocate(tmpVec(nDim))
    tmpVec = 0
    allocate(tmpMat(nDim,nNucOrb))
    tmpMat = 0
 
    do i = 1,nDim
      nLapl(i,i) = -73766.0/25200.0
      if (i < nDim) then 
        nLapl(i,i+1) = 5.0/3.0
      endif
      if (i > 1) then
          nLapl(i,i-1) = ( 5.0/3.0)
      endif
      if (i < (nDim - 1)) then
          nLapl(i,i+2) = ( -5.0/21.0)
      endif
      if (i > 2) then
          nLapl(i,i-2) = ( -5.0/21.0)
      endif
      if (i < (nDim - 2)) then
          nLapl(i,i+3) = ( 5.0/126.0)
      endif
      if (i > 3) then
          nLapl(i,i-3) = ( 5.0/126.0)
      endif
    if (i < (nDim - 3)) then
          nLapl(i,i+4) = ( -5.0/1008.0)
      endif
      if (i > 4) then
          nLapl(i,i-4) = ( -5.0/1008.0)
      endif
      if (i < (nDim - 4)) then
          nLapl(i,i+5) = ( +1.0/3150.0)
      endif
      if (i > 5) then
          nLapl(i,i-5) = ( +1.0/3150.0)
      endif    
    enddo
 
    Tn = -(1/(2*dxn**2*mu_n))*nLapl
  
    !!!!!!!!!!!!!!!!!!!!! Nuclear Basis   !!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(chi(nDim,nNucOrb))
    chi=0
    allocate(eigenVals(nNucOrb))
    eigenVals=0
    call diagSymMatrix(Tn, nNucOrb, eigenVals, chi) 
    
    Tn = -(1/(2*dxn**2*mu_n))*nLapl

    do i=1,nNucOrb
      chi(:,i) = chi(:,i)/sqrt(sum(chi(:,i)*chi(:,i))*dxn)
    enddo

    !!!!!!!!!!!!!!!!!!!!! Nuclear Kernel Matrices !!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! H = Tn + Te + Vee + Vnn + Ven

    allocate(KTn(nNucOrb,nNucOrb))
    KTn = 0
    allocate(Vnn(nDim,nDim))
    Vnn = 0
    allocate(KVnn(nNucOrb,nNucOrb))
    KVnn = 0
    do i = 1, nDim 
      Vnn(i,i) = 1.0/nAxis(i)
    end do
    !Calculate Vnn Matrix 
    call sgemm('N', 'N', nDim, nNucOrb, nDim, 1.0, Vnn, nDim, chi, nDim, 0, tmpMat, nDim)
    
    !tmpMat is now Vnn*chi
    call sgemm('T', 'N', nNucOrb, nNucOrb, nDim, dxn, chi, nDim, tmpMat, nDim, 0, KVnn, nNucOrb)
    
    !tmpMat = 0
    !Calculate VTn Matrix 
    call sgemm('N', 'N', nDim, nNucOrb, nDim, 1.0, Tn, nDim, chi, nDim, 0, tmpMat, nDim)
    
    !tmpMat is now Tn*chi
    call sgemm('T', 'N', nNucOrb, nNucOrb, nDim, dxn, chi, nDim, tmpMat, nDim, 0, KTn, nNucOrb)

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Load Electronic Information !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !TODO: Replace all of these with sparse reads
 
    ! This matrix is printed in lower diagonal form, simply write as upper diag 
    TeFile = 1
    allocate(KTe(NeleSpatialOrb,NeleSpatialOrb))
    KTe = 0
    open(TeFile, file='RL_0.50_RR_6.00_dR_0.10_NEx_1/Te_Vee/Te', status='old')
    
    do while (continueReading)
      read(TeFile,*) i, j, zval
     
      if(imagpart(zval) > 0 ) then
        error stop
      endif
       
      ! This works because i increases slowest
      if ( i <= NeleSpatialOrb .and. j <= NeleSpatialOrb ) then
        KTe(j,i) = realpart(zval) !note j, i to make upper triangular
        if ( i==NeleSpatialOrb .and. j==NeleSpatialOrb ) then
          continueReading = .false.
        endif
      endif
    end do

    close(TeFile)
    !print *, "Te is: "
    !do i=1,NeleSpatialOrb
    !  print *, KTe(i,:)
    !enddo

    VeeFile = 2
    allocate(KVee(NeleSpatialOrb**2,NeleSpatialOrb**2))
    KVee = 0
    open(VeeFile,file='RL_0.50_RR_6.00_dR_0.10_NEx_1/Te_Vee/Vee', status='old')
    ! Dummy read of header line
    read(VeeFile,*) 

    continueReading = .true.
    do while (continueReading)
      !              n1,k1,n2,k2,n3,k3,n4,k4, val 
      read(VeeFile,*) i, j, k, l, m, n, o, p, rval
     
      
      ! If we've read in all of the spatial orbitals needed, stop 
      if ( i <= NeleSpatialOrb .and. k <= NeleSpatialOrb  .and. m <= NeleSpatialOrb .and. o<= NeleSpatialOrb ) then
        ! Fold the Vee spatial integrals 
        ii = k + NeleSpatialOrb*(i-1) 
        jj = o + NeleSpatialOrb*(m-1)
        KVee(ii,jj) = rval
        if ( i==NeleSpatialOrb .and. k==NeleSpatialOrb .and. m==NeleSpatialOrb .and. o==NeleSpatialOrb ) then
          continueReading = .false.
        endif
      endif
    end do
    close(VeeFile)

    VenFile = 3
    allocate(KVen(nDim, NeleSpatialOrb, NeleSpatialOrb))
    KVen = 0


    !This matrix is written as lower diagonal, simply write as upper diag
    !TODO Adjust for R values over 9.99
    do Ri=1,nDim
      write(filepath, "(f4.2)") nAxis(Ri)
      filepath = "RL_0.50_RR_6.00_dR_0.10_NEx_1/Ven/R_"//filepath
     
      filepath = trim(filepath)//"/output_me_one_body"

      continueReading=.true.
      open(VenFile,file=filepath,status='old')
      
      do while (continueReading)
              
        read(VenFile,*) i, j, zval
        
        if(imagpart(zval) > 0 ) then
          error stop
        endif
         
        ! This works because i increases slowest
        if ( i <= NeleSpatialOrb .and. j <= NeleSpatialOrb ) then
          KVen(Ri,j,i) = realpart(zval) !note j, i to make upper triangular
          if ( i==NeleSpatialOrb .and. j==NeleSpatialOrb ) then
            continueReading = .false.
          endif
        endif 
      enddo !read
    enddo ! R scan

    ! The matrix elements are arranged with the fastest index being the nuclear orbitals
    ! and the slowest index being the slater determinants

    allocate(spinIndI(Nele))
    spinIndI = 0
    allocate(spaceIndI(Nele))
    spaceIndI = 0
    allocate(spinIndIp(Nele))
    spinIndI = 0
    allocate(spaceIndIp(Nele))
    spaceIndI = 0

    allocate(nucRange(NnucOrb))
    nucRange = (/ (i, i=1,NnucOrb) /)
    allocate(slaterIJ(NnucOrb))
    allocate(slaterIpJp(NnucOrb))
    allocate(MTn(NSlater*NnucOrb,NSlater*NnucOrb))
    allocate(MVee(NSlater*NnucOrb,NSlater*NnucOrb))
    allocate(MTe(NSlater*NnucOrb,NSlater*NnucOrb))
    allocate(IdN(NnucOrb,NnucOrb)) 
    IdN = 0
    do i=1,NnucOrb
      IdN(i,i) = 1.0
    enddo

    do slaterI=1,NSlater
      do slaterIp=slaterI,NSlater
        slaterIJ = nucRange + NnucOrb*(slaterI-1)
        slaterIpJp = nucRange + NnucOrb*(slaterIp-1)

        !MTn
        if(slaterIp==slaterI) then !TODO: only upper triangular assign
          MTn(slaterIJ,slaterIpJp) = KTn
        endif

        
        if(CITruncation=='doubles') then
          !In this case there are no slater determinants which vary by
          !only one spin orbital
            
          !Get spin orbital indeces associated to Ip. 
          spinIndI    =   slaterIndices(slaterI)%spinOrbitals;
          spaceIndI   =   slaterIndices(slaterI)%spatialOrbitals;
          spinIndIp   =   slaterIndices(slaterIp)%spinOrbitals;
          spaceIndIp  =   slaterIndices(slaterIp)%spatialOrbitals;
          if(slaterIp==slaterI) then
            !TODO: 14th Sep check this
            !MTe 
            MTeVal=0;
            MVeeVal =0;
            do i=1,Nele
                ! Szabo Table 2.3 Case 1
                MTeVal = MTeVal + KTe(spaceIndI(i),spaceIndI(i));
                do j=(i+1),Nele 
                  ! j+1 as the exchange makes identical spatial orbital integrals die
                  ! in this case
                  ii = spaceIndI(i) + NeleSpatialOrb*(spaceIndI(i)-1);
                  jj = spaceIndI(j) + NeleSpatialOrb*(spaceIndI(j)-1);
                  ij = spaceIndI(i) + NeleSpatialOrb*(spaceIndI(j)-1);
                  ji = spaceIndI(j) + NeleSpatialOrb*(spaceIndI(i)-1);
            
                  ! Szabo Table 2.4 Case 1
                  if(mod(spinIndI(i),2) == mod(spinIndI(j),2)) then
                    MVeeVal = MVeeVal  +  KVee(ii,jj) - KVee(ij,ji);
                  else
                    MVeeVal = MVeeVal  +  KVee(ii,jj);
                  endif
                enddo
            enddo
            MTe(slaterIJ,slaterIpJp) = MTeVal*IdN; 
            MVee(slaterIJ,slaterIpJp) = MVeeVal*IdN;
          else !if not on diagonal 
            ! We have to put these into maximum coincidence
 
          endif ! slaterI==slaterIP, else
        endif ! CITruncation=='doubles'
      enddo ! SlaterIp 
    enddo ! Slater I
    
    print *, "MVee is : "
    do i=1,NnucOrb*NSlater
      print *, MVee(i,:)
    enddo

    
  end subroutine setOperators   


end module operators
