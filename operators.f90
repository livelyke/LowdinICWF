module operators
  use params, only          : mu_n, dxn, nDim, Nele, NeleSpinOrb, NSlater, NeleSpatialOrb, &
                              NnucOrb, nAxis, CITruncation, slaterIndices, debug, electronicOnly
  use bases, only           : diagSymMatrix
  use modSlaterIndex, only  : SlaterIndex, SlaterMaximumCoincidence

  implicit none (type, external)
  external  :: ssyev, sgemm, sgeev, ssyevr, amux
  private
  public    :: setOperators

  real, dimension(:,:), allocatable, public   ::  H

  contains

  ! TODO: Naming convention here is a bit shit
  !       replace spaceInd / spinInd with spaceOrb
  !       change name of dSpaceOrbI1 potentially
  subroutine setOperators()
    integer, dimension(:), allocatable          :: nucRange, slaterIJ, slaterIpJp
    real, dimension(:,:), allocatable           ::  chi, Tn, nLapl, Vnn, KTn, KVnn, &
                                                    IdN, IdE, KTe, KVee, MTn, MTe, MVee, &
                                                    MVen, MVnn, tmpMat1, tmpMat2, tmpMat, MVenVals, &
                                                    spinOrbsI, spinOrbsIp, spaceOrbsI, spaceOrbsIp
    real, dimension(:,:,:), allocatable         ::  KVen
    real, dimension(:), allocatable             ::  eigenVals, VenKernalOp, coefsI, coefsIp
    real                                ::  MTeVal=0, MVeeVal=0, tmpVal1, tmpVal2
    integer                             ::  i, j, m, n, &
                                            slaterI, slaterIp, &
                                            sgn, indexLocation, &
                                            numDiffering, dSpaceOrbI1, dSpaceOrbI2, dSpaceOrbIp1, &
                                            dSpaceOrbIp2, &
                                            im_im, in_in, im_in, in_im, im_ip, in_ip, in_iq, im_iq, &
                                            dSpinOrbI1, dSpinOrbIp1, dSpinOrbI2, dSpinOrbIp2, &
                                            NConfI, NConfIp


    !TODO: check if rearranged nuclear set up here still agrees with matlab
    !>>>>> Set Nuclear Operators <<<<<!
    allocate(Tn(nDim,nDim))
    Tn=0
    allocate(nLapl(nDim,nDim))
    nLapl=0
    call setNuclearOperators(Tn,nLapl)

    !>>>>> Set Nuclear Basis <<<<<!
    allocate(eigenVals(nNucOrb))
    eigenVals=0
    allocate(chi(nDim,nNucOrb))
    chi = 0
    call diagSymMatrix(Tn, nNucOrb, eigenVals, chi) 
    
    Tn = -(1/(2*dxn**2*mu_n))*nLapl

    do i=1,nNucOrb
      chi(:,i) = chi(:,i)/sqrt(sum(chi(:,i)*chi(:,i))*dxn)
    enddo

    !>>>>> Nuclear Kernel Matrices <<<<<!
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

    allocate(tmpMat(nDim,NnucOrb))
    tmpMat = 0
    !Calculate Vnn Matrix 
    call sgemm('N', 'N', nDim, nNucOrb, nDim, 1.0, Vnn, nDim, chi, nDim, 0, tmpMat, nDim)
    
    !tmpMat is now Vnn*chi
    call sgemm('T', 'N', nNucOrb, nNucOrb, nDim, dxn, chi, nDim, tmpMat, nDim, 0, KVnn, nNucOrb)
    
    !tmpMat = 0
    !Calculate VTn Matrix 
    call sgemm('N', 'N', nDim, nNucOrb, nDim, 1.0, Tn, nDim, chi, nDim, 0, tmpMat, nDim)
    
    !tmpMat is now Tn*chi
    call sgemm('T', 'N', nNucOrb, nNucOrb, nDim, dxn, chi, nDim, tmpMat, nDim, 0, KTn, nNucOrb)

    
    !>>>>>>>>>>>>>>> Load Electronic Information <<<<<<<<<<<<<<<<<!
 
    allocate(KTe(NeleSpatialOrb,NeleSpatialOrb))
    KTe = 0
    call readTe('Te',KTe)

    allocate(KVee(NeleSpatialOrb**2,NeleSpatialOrb**2))
    KVee = 0
    call readVee('Vee',KVee)

    allocate(KVen(nDim, NeleSpatialOrb, NeleSpatialOrb))
    KVen = 0
    call readVen('Ven',KVen)
   
    !>>>>>>>>>>>>>>> Set Hamiltonian Matrix Elements <<<<<<<<<<<<<<<! 
    ! The matrix elements are arranged with the fastest index being the nuclear orbitals
    ! and the slowest index being the slater determinants

    allocate(spinOrbsI(Nele, 6))
    spinOrbsI = 0
    allocate(spaceOrbsI(Nele, 6))
    spaceOrbsI = 0
    allocate(spinOrbsIp(Nele, 6))
    spinOrbsIp = 0
    allocate(spaceOrbsIp(Nele, 6))
    spaceOrbsIp = 0

    allocate(coefsI(6))
    coefsI = 0
    allocate(coefsIp(6))
    coefsIp = 0

    allocate(nucRange(NnucOrb))
    nucRange = (/ (i, i=1,NnucOrb) /)
    allocate(slaterIJ(NnucOrb))
    allocate(slaterIpJp(NnucOrb))
    allocate(MTn(NSlater*NnucOrb,NSlater*NnucOrb))
    MTn = 0
    allocate(MVnn(NSlater*NnucOrb,NSlater*NnucOrb))
    MVnn = 0
    allocate(MVee(NSlater*NnucOrb,NSlater*NnucOrb))
    MVee = 0
    allocate(MTe(NSlater*NnucOrb,NSlater*NnucOrb))
    MTe = 0
    allocate(MVen(NSlater*NnucOrb,NSlater*NnucOrb))
    MVen = 0
    allocate(MVenVals(NnucOrb,NnucOrb))
    MVenVals = 0
    allocate(VenKernalOp(nDim))
    VenKernalOp = 0
    allocate(IdN(NnucOrb,NnucOrb)) 
    IdN = 0
    do i=1,NnucOrb
      IdN(i,i) = 1.0
    enddo
    allocate(IdE(NSlater,NSlater)) 
    IdE = 0
    do i=1,NSlater
      IdE(i,i) = 1.0
    enddo

    if (debug .eqv. .true.) then
      allocate(tmpMat1(NSlater,NSlater))
      allocate(tmpMat2(NSlater,NSlater))
       tmpMat1 = 0
       tmpMat2 = 0
    endif

    !TODO: double check that MVen is working for non-electronicOnly
    do slaterI=1,NSlater
      do slaterIp=slaterI,NSlater
        
        slaterIJ = nucRange + NnucOrb*(slaterI-1)
        slaterIpJp = nucRange + NnucOrb*(slaterIp-1)

        !MTn
        if(slaterIp==slaterI) then !TODO: only upper triangular assign, in sparse fashion
          MTn(slaterIJ,slaterIpJp) = KTn
          MVnn(slaterIJ,slaterIpJp) = KVnn
        endif
        
        !> Get spin orbital indices associated to Ip. 
        NConfI  = slaterIndices(slaterI)%NConf
        NConfIp = slaterIndices(slaterI)%NConf

        spinOrbsI(:,1:NConfI)   = slaterIndices(slaterI)%spinOrbitals;
        spaceOrbsI(:,1:NConfI)  = slaterIndices(slaterI)%spaceOrbitals;
        spinOrbsIp(:,1:NConfIp)  = slaterIndices(slaterIp)%spinOrbitals;
        spaceOrbsIp(:,1:NConfIp) = slaterIndices(slaterIp)%spaceOrbitals;

        MTeVal  = 0;
        MVeeVal = 0;
        MVenVals = 0;
        VenKernalOp = 0;
        
        do i=1,NConfI
          do j=1,NConfIp
          enddo
        enddo

      enddo ! SlaterIp 
    enddo ! Slater I

    if (debug .eqv. .true.) then
      print *, "MTn"
      do i=1,NSlater*NnucOrb
        print *, MTn(i,:)
      enddo
      print *, "MTe"
      do i=1,NSlater*NnucOrb
        print *, MTe(i,:)
      enddo
      print *, "MVee Unexchanged"
      do i=1,NSlater*NnucOrb
        print *, tmpMat1(i,:)
      enddo
      print *, "MVee Exchanged"
      do i=1,NSlater*NnucOrb
        print *, tmpMat2(i,:)
      enddo
      print *, "MVee"
      do i=1,NSlater*NnucOrb
        print *, MVee(i,:)
      enddo
      print *, "MVen"
      do i=1,NSlater*NnucOrb
        print *, MVen(i,:)
      enddo
      print *, "MVnn"
      do i=1,NSlater*NnucOrb
        print *, MVnn(i,:)
      enddo
    endif

    if(electronicOnly) then
      allocate(H(NSlater,NSlater))
      H = 0
      H = MTe + MVee + MVen + (Vnn(1,1))*IdE
    else
      H = MTn + MTe + MVee + MVen + MVnn
    endif
   
    deallocate(Tn) 
    deallocate(nLapl)
    deallocate(eigenVals)
    deallocate(chi)
    deallocate(KTn) 
    deallocate(Vnn) 
    deallocate(KVnn) 
    deallocate(KVee) 
    deallocate(KVen) 
    deallocate(KTe)
    deallocate(spinOrbsI) 
    deallocate(spaceOrbsI) 
    deallocate(spinOrbsIp) 
    deallocate(spaceOrbsIp)
    deallocate(nucRange)
    deallocate(slaterIJ) 
    deallocate(slaterIpJp)
    deallocate(MTn) 
    deallocate(MTe) 
    deallocate(MVee) 
    deallocate(MVen) 
    deallocate(MVenVals) 
    deallocate(VenKernalOp) 
    deallocate(IdN) 
    deallocate(IdE) 
    deallocate(MVnn) 
    deallocate(tmpMat)
    if(debug .eqv. .true.) then  
      deallocate(tmpMat1) 
      deallocate(tmpMat2) 
    endif
  end subroutine setOperators   

  subroutine setNuclearOperators(Tn, nLapl)
    real, intent(inout)                 ::  Tn(:,:), nLapl(:,:)
    integer         ::  i

    nLapl = 0
    Tn = 0
 
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
  
  end subroutine setNuclearOperators

  !TODO: Replace all of these with sparse reads
  subroutine readTe(filePath, KTe)
    character(len=*), intent(in)     :: filePath
    real,               intent(out)    :: KTe(:,:)
    integer                            :: TeFile,i,j
    complex                            :: zval
    logical                            :: continueReading=.true. 
    
    ! This matrix is printed in lower diagonal form, simply write as upper diag 
    TeFile = 1
    open(TeFile, file=trim(filePath), status='old')
    
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
    
    if (debug .eqv. .true.) then
      print *, "Te is: "
      do i=1,NeleSpatialOrb
        print *, KTe(i,:)
      enddo
    endif

  end subroutine readTe 

  subroutine readVee(filePath,KVee)
    character(len=*), intent(in)     :: filePath
    real,               intent(out)    :: KVee(:,:)
    integer                            :: VeeFile,i,j,k,l,m,n,o,p,ii,jj
    real                               :: rval
    logical                            :: continueReading=.true. 
    
    VeeFile = 2
    open(VeeFile,file=trim(filePath), status='old')
    ! Dummy read of header line
    read(VeeFile,*) 

    continueReading = .true.
    do while (continueReading)
      !              n1,k1,n2,k2,n3,k3,n4,k4, val 
      read(VeeFile,*) i, j, k, l, m, n, o, p, rval
     
      
      ! If we've read in all of the spatial orbitals needed, stop 
      if ( i <= NeleSpatialOrb .and. k <= NeleSpatialOrb  .and. m <= NeleSpatialOrb .and. o<= NeleSpatialOrb ) then
        ! Fold the Vee spatial integrals
        ! Fold such that for (n1 n2 | n3 n4)
        ! n1 is fastest and n3 is fastest 
        ii = i + NeleSpatialOrb*(k-1) 
        jj = m + NeleSpatialOrb*(o-1)
        KVee(ii,jj) = rval
        if ( i==NeleSpatialOrb .and. k==NeleSpatialOrb .and. m==NeleSpatialOrb .and. o==NeleSpatialOrb ) then
          continueReading = .false.
        endif
      endif
    end do
    close(VeeFile)
    
    if (debug .eqv. .true.) then
      print *, "Vee is: "
      do i=1,NeleSpatialOrb
        print *, KVee(i,:)
      enddo
    endif

  end subroutine readVee

  subroutine readVen(filePathBase, KVen)
    character(len=*), intent(in)        :: filePathBase
    real,               intent(out)     :: KVen(:,:,:)
    character(len=256)                    :: rstring, filePath
    complex                             :: zval
    logical                             :: continueReading=.true.
    integer                             :: i,j,VenFile, Ri 

    VenFile = 3
    !This matrix is written as lower diagonal, simply write as upper diag
    !TODO: Adjust for R values over 9.99
    !TODO: double check
    do Ri=1,nDim
      write(rstring, "(f4.2)") nAxis(Ri)
      rstring = trim(rstring)
      filePath = trim(filePathBase)//"/R_"
      filePath = trim(filePath)//rstring
      filePath = trim(filePath)//"/output_me_one_body"
    
      continueReading=.true.
      open(VenFile,file=filePath,status='old')
      
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

    if (debug .eqv. .true.) then
      open(1,file='Ven.debug',status='old')
      do i=1,nDim
          write(1,*) nAxis(i)
          do j=1,NeleSpatialOrb
            write(1,*) KVen(i,j,:)
          enddo
      enddo
      close(1)
    endif
  end subroutine readVen

end module operators
