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
    integer, dimension(:), allocatable    :: nucRange, slaterIJ, slaterIpJp
    integer, dimension(:,:), allocatable  :: spinOrbsI, spinOrbsIp, spaceOrbsI, spaceOrbsIp
    real, dimension(:,:), allocatable     :: chi, Tn, nLapl, Vnn, KTn, KVnn, &
                                             IdN, IdE, KTe, KVee, MTn, MTe, MVee, &
                                             MVen, MVnn, tmpMat1, tmpMat2, tmpMat, MVenVals, &
                                             MVenValsTmp, KDx, KDy, KDz
    real, dimension(:,:,:), allocatable   :: KVen
    real, dimension(:), allocatable       :: eigenVals, VenKernalOp, VenKernalOpTmp, coefsI, coefsIp
    real                                  :: MTeVal, MVeeVal, MTeValTmp, MVeeValTmp
    integer                               :: i, j, m, n, &
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
    
    allocate(KDx(NeleSpatialOrb,NeleSpatialOrb))
    KDx = 0
    allocate(KDy(NeleSpatialOrb,NeleSpatialOrb))
    KDy = 0
    allocate(KDz(NeleSpatialOrb,NeleSpatialOrb))
    KDz = 0

    call readDipole('ks_me_dipole.k1_',KDx,KDy,KDz)   

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
    allocate(MVenValsTmp(NnucOrb,NnucOrb))
    MVenValsTmp = 0
    allocate(VenKernalOp(nDim))
    VenKernalOp = 0
    allocate(VenKernalOpTmp(nDim))
    VenKernalOpTmp = 0
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
        NConfIp = slaterIndices(slaterIp)%NConf

        spinOrbsI(:,1:NConfI)    = slaterIndices(slaterI)%spinOrbitals;
        spaceOrbsI(:,1:NConfI)   = slaterIndices(slaterI)%spaceOrbitals;
        coefsI(1:NConfI)         = slaterIndices(slaterI)%coefs;
        spinOrbsIp(:,1:NConfIp)  = slaterIndices(slaterIp)%spinOrbitals;
        spaceOrbsIp(:,1:NConfIp) = slaterIndices(slaterIp)%spaceOrbitals;
        coefsIp(1:NConfIp)       = slaterIndices(slaterIp)%coefs;

        MTeVal      = 0
        MVeeVal     = 0
        MVenVals    = 0
        VenKernalOp = 0
        
        do i=1,NConfI
          do j=1,NConfIp
            if (debug .eqv. .true.) then
              write(*,'(A,I3,A,I1,A,I3,A,I1)') "I=", slaterI, ".",i,"  Ip=",slaterIp,".",j
            endif
            !> Put the constituent slater determinants into maximum coincidence
            call SlaterMaximumCoincidence(          &
                  spinOrbsI(:,i),spinOrbsIp(:,j),   &
                  sgn, numDiffering, dSpinOrbI1,     &
                  dSpinOrbI2, dSpinOrbIp1, dSpinOrbIp2)

            if (numDiffering <= 1) then
              MTeValTmp      = 0
              MVeeValTmp     = 0
              MVenValsTmp    = 0
              VenKernalOpTmp = 0
            elseif (numDiffering == 2) then
              MVeeValTmp     = 0
            endif

            if (numDiffering == 0) then
              call matrixElementsNull(        &
                                MTeValTmp, MVeeValTmp, MVenValsTmp, VenKernalOpTmp, &
                                KTe, KVen, KVee,&
                                spaceOrbsI(:,i), spinOrbsI(:,i))
            elseif (numDiffering == 1) then
              if (debug .eqv. .true.) then
                write(*,'(A,I3,A,I1,A,I3,A,I1)') "numDiffering 1, I=", slaterI, ".",i,"  Ip=",slaterIp,".",j
              endif
              call matrixElementsOne(        &
                                MTeValTmp, MVeeValTmp, MVenValsTmp, VenKernalOpTmp, &
                                KTe, KVen, KVee,&
                                spaceOrbsI(:,i), spinOrbsI(:,i), spaceOrbsIp(:,j),spinOrbsIp(:,j), &
                                dSpinOrbI1, dSpinOrbIp1)
            elseif (numDiffering == 2) then
              if (debug .eqv. .true.) then
                write(*,'(A,I3,A,I1,A,I3,A,I1)') "numDiffering 2, I=", slaterI, ".",i,"  Ip=",slaterIp,".",j
                    print *, "<", dSpinOrbI1, dSpinOrbI2, "  | |",dSpinOrbIp1, dSpinOrbIp2, ">"
              endif
              call matrixElementsTwo(        &
                                MVeeValTmp, KVee,&
                                spaceOrbsI(:,i), spinOrbsI(:,i), spaceOrbsIp(:,j),spinOrbsIp(:,j), &
                                dSpinOrbI1, dSpinOrbIp1, dSpinOrbI2, dSpinOrbIp2)
            endif
            !> Sign from Maximum Coincidence is applied here
            if (numDiffering <= 1) then
              !> If numDiffering = 0 then MTeValTmp,MVenValsTmp, and VenKernalOpTmp = 0
              MTeVal      = MTeVal      + sgn*coefsI(i)*coefsIp(j)*MTeValTmp
              MVeeVal     = MVeeVal     + sgn*coefsI(i)*coefsIp(j)*MVeeValTmp
              MVenVals    = MVenVals    + sgn*coefsI(i)*coefsIp(j)*MVenValsTmp
              VenKernalOp = VenKernalOp + sgn*coefsI(i)*coefsIp(j)*VenKernalOpTmp
            elseif (numDiffering == 2) then
              MVeeVal     = MVeeVal     + sgn*coefsI(i)*coefsIp(j)*MVeeValTmp
            endif
          enddo !< configurations of slater Ip
        enddo !< configurations of slater I

        MTe(slaterIJ,slaterIpJp) = MTeVal*IdN; 
        MVee(slaterIJ,slaterIpJp) = MVeeVal*IdN;
        if (.not. electronicOnly) then
          do i=1,NnucOrb
            do j=i,NnucOrb
              MVenVals(i,j) = sum(chi(:,i)*VenKernalOp*chi(:,j))*dxn
            enddo
          enddo
        endif
        MVen(slaterIJ,slaterIpJp) = MVenVals
        
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


  subroutine matrixElementsNull(MTeVal, MVeeVal, MVenVals, VenKernalOp, &
                                KTe, KVen, KVee,&
                                spaceOrbsI, spinOrbsI)

    real, intent(inout)           :: MTeVal, MVeeVal, MVenVals(:,:), VenKernalOp(:), &
                                     KTe(:,:), KVen(:,:,:), KVee(:,:)
    integer, intent(in)           :: spaceOrbsI(:), spinOrbsI(:)
    integer   ::  m, n, im_im, in_in, im_in, in_im

    
    do m=1,Nele
        ! Szabo Table 2.3 Case 1
        MTeVal = MTeVal + KTe(spaceOrbsI(m),spaceOrbsI(m));
        if (electronicOnly) then
          MVenVals = MVenVals + KVen(1,spaceOrbsI(m),spaceOrbsI(m))
        else
          VenKernalOp = VenKernalOp + KVen(:,spaceOrbsI(m),spaceOrbsI(m))
        endif

        do n=(m+1),Nele 
          ! j+1 as the exchange makes identical spatial orbital integrals die
          ! for slaterI == slaterIp
          im_im = spaceOrbsI(m) + NeleSpatialOrb*(spaceOrbsI(m)-1);
          in_in = spaceOrbsI(n) + NeleSpatialOrb*(spaceOrbsI(n)-1);
          im_in = spaceOrbsI(m) + NeleSpatialOrb*(spaceOrbsI(n)-1);
          in_im = spaceOrbsI(n) + NeleSpatialOrb*(spaceOrbsI(m)-1);
    
          ! Szabo Table 2.4 Case 1
          if(mod(spinOrbsI(m),2) == mod(spinOrbsI(n),2)) then
            MVeeVal = MVeeVal  +  KVee(im_im,in_in) - KVee(im_in,in_im);
          else
            MVeeVal = MVeeVal  +  KVee(im_im,in_in);
          endif
        enddo
    enddo

  end subroutine matrixElementsNull

  subroutine matrixElementsOne(MTeVal, MVeeVal, MVenVals, VenKernalOp, &
                               KTe, KVen, KVee,&
                               spaceOrbsI, spinOrbsI, spaceOrbsIp, spinOrbsIp, &
                               dSpinOrbI1, dSpinOrbIp1)

    real, intent(inout)           :: MTeVal, MVeeVal, MVenVals(:,:), VenKernalOp(:), &
                                     KTe(:,:), KVen(:,:,:), KVee(:,:)
    integer, intent(in)           :: spaceOrbsI(:), spinOrbsI(:), spinOrbsIp(:), spaceOrbsIp(:), &
                                     dSpinOrbI1, dSpinOrbIp1
    integer   ::  m, n, im_ip, in_in, im_in, in_ip, dSpaceOrbI1, dSpaceOrbIp1, indexLocation

 
    indexLocation = findloc(spinOrbsI,dSpinOrbI1,1)
    dSpaceOrbI1    = spaceOrbsI(indexLocation)
    
    indexLocation = findloc(spinOrbsIp,dSpinOrbIp1,1)
    dSpaceOrbIp1   = spaceOrbsIp(indexLocation)
 
    ! Szabo Table 2.3 Case 2
    ! KTe and KVen are symmetric and stored as upper triangular
    if( dSpaceOrbI1 > dSpaceOrbIp1 ) then
      MTeVal = KTe(dSpaceOrbIp1,dSpaceOrbI1)
      VenKernalOp = KVen(:,dSpaceOrbIp1, dSpaceOrbI1) 
    else
      MTeVal = KTe(dSpaceOrbI1,dSpaceOrbIp1)
      VenKernalOp = KVen(:,dSpaceOrbI1, dSpaceOrbIp1) 
    endif
 
    ! Szabo Table 2.4 Case 2
    ! see also 2.148
    ! sum_{m\neq n } <mn| |pn> 
    ! sum_{m\neq n} (i_m i_p | i_n i_n)\delta_{m%2, p%2} - (i_m i_n | i_n i_p)\delta_{m%2,n%2}\delta_{n%2,p%2} 
    ! recall that KVee is folded for (n1 n2 | n3 n4)
    ! as i_n1 + NeleSpatialOrb*(i_n2 - 1)
    !    i_n3 + NeleSpatialOrb*(i_n4 - 1)
    !if (debug .eqv. .true.) then
    !  tmpVal1 = 0
    !  tmpVal2 = 0
    !endif
    
    ! in the naming convention so far we call differing m, p as dSpaceOrbI1, dSpaceOrbIp1
    im_ip = dSpaceOrbI1 + NeleSpatialOrb*(dSpaceOrbIp1-1)
    do n=1,Nele
      ! \sum_n <mn| |pn> 
      ! except <mm| |pm> which is zero
      if ( spinOrbsI(n) .ne. dSpinOrbI1) then  
        if (debug .eqv. .true.) then
          print *, "<",dSpinOrbI1, spinOrbsI(n), "  | |",dSpinOrbIp1, spinOrbsI(n), ">"
        endif
    
        in_in = spaceOrbsI(n) + NeleSpatialOrb*(spaceOrbsI(n)-1)
        if ( mod(dSpinOrbI1,2) .eq. mod(dSpinOrbIp1,2) ) then
          MVeeVal = MVeeVal + KVee(im_ip,in_in)
          if (debug .eqv. .true.) then
            print *, "(", dSpaceOrbI1, dSpaceOrbIp1, "   | ", spaceOrbsI(n), spaceOrbsI(n), ") = ", KVee(im_ip,in_in)
            !tmpVal1 = tmpVal1 + KVee(im_ip,in_in)
          endif
        endif
    
        if ( (mod(dSpinOrbI1,2) .eq. mod(spinOrbsI(n),2))  .and. (mod(dSpinOrbIp1,2) .eq. mod(spinOrbsI(n),2)) ) then
          im_in = dSpaceOrbI1   + NeleSpatialOrb*(spaceOrbsI(n)- 1)
          in_ip = spaceOrbsI(n) + NeleSpatialOrb*(dSpaceOrbIp1 - 1)
          MVeeVal = MVeeVal - KVee(im_in,in_ip)
          if (debug .eqv. .true.) then
            print *, "(", dSpaceOrbI1, spaceOrbsI(n), "   | ", spaceOrbsI(n), dSpaceOrbIp1, ") = ", KVee(im_in,in_ip)
            !tmpVal2 = tmpVal2 - KVee(im_in,in_ip)
          endif
        endif
    
        !>>> Add the appropriate sign <<<!
        !if (debug .eqv. .true.) then
        !  tmpVal1 = sgn*tmpVal1
        !  tmpVal2 = sgn*tmpVal2
        !endif
    
      endif ! check if <mm| |pm> = 0
    enddo ! Loop over spin Orbitals
    
  end subroutine matrixElementsOne

  subroutine matrixElementsTwo(MVeeVal, KVee,&
                               spaceOrbsI, spinOrbsI, spaceOrbsIp, spinOrbsIp, &
                               dSpinOrbI1, dSpinOrbIp1, dSpinOrbI2, dSpinOrbIp2)

    real, intent(inout)           :: MVeeVal, KVee(:,:)
    integer, intent(in)           :: spaceOrbsI(:), spinOrbsI(:), spinOrbsIp(:), spaceOrbsIp(:), &
                                     dSpinOrbI1, dSpinOrbIp1, dSpinOrbI2, dSpinOrbIp2
    integer   ::  im_ip, in_iq, im_iq, in_ip, dSpaceOrbI1, dSpaceOrbIp1, &
                  dSpaceOrbI2, dSpaceOrbIp2, indexLocation

    ! <mn| |pq> = <mn||pq> - <mn|qp>
    ! m = dSpinOrbI1
    ! n = dSpinOrbI2
    ! p = dSpinOrbIp1
    ! q = dSpinOrbIp2
    
    indexLocation = findloc(spinOrbsI,dSpinOrbI1,1)
    dSpaceOrbI1    = spaceOrbsI(indexLocation)
    indexLocation = findloc(spinOrbsI,dSpinOrbI2,1)
    dSpaceOrbI2    = spaceOrbsI(indexLocation)
    
    indexLocation = findloc(spinOrbsIp,dSpinOrbIp1,1)
    dSpaceOrbIp1   = spaceOrbsIp(indexLocation)
    indexLocation = findloc(spinOrbsIp,dSpinOrbIp2,1)
    dSpaceOrbIp2   = spaceOrbsIp(indexLocation)
    
    !if (debug .eqv. .true.) then
    !  tmpVal1 = 0
    !  tmpVal2 = 0
    !endif
    if( (mod(dSpinOrbI1,2) .eq. mod(dSpinOrbIp1,2)) .and. (mod(dSpinOrbI2,2) .eq. mod(dSpinOrbIp2,2)) ) then
    
      ! <mn||pq> = (im_ip|in_iq) \delta_{m%2,p%2}\delta_{n%2,q%2}
      
      im_ip = dSpaceOrbI1 + NeleSpatialOrb*(dSpaceOrbIp1 - 1)
      in_iq = dSpaceOrbI2 + NeleSpatialOrb*(dSpaceOrbIp2 - 1)
      
      MVeeVal = MVeeVal + KVee(im_ip, in_iq)
      if (debug .eqv. .true.) then
        print *, "(", dSpaceOrbI1, dSpaceOrbIp1, "   | ",dSpaceOrbI2, dSpaceOrbIp2, ") = ", KVee(im_ip,in_iq)
        !tmpVal1 = tmpVal1 + KVee(im_ip, in_iq)
      endif
    
    endif
    if( (mod(dSpinOrbI1,2) .eq. mod(dSpinOrbIp2,2)) .and. (mod(dSpinOrbI2,2) .eq. mod(dSpinOrbIp1,2)) ) then
      
      ! -1*<mn||qp> = -1*(im_iq|in_ip)\delta_{m%2,q%2}\delta_{n%2,p%2}
      
      im_iq = dSpaceOrbI1 + NeleSpatialOrb*(dSpaceOrbIp2 - 1)
      in_ip = dSpaceOrbI2 + NeleSpatialOrb*(dSpaceOrbIp1 - 1)
      
      MVeeVal = MVeeVal - KVee(im_iq, in_ip)
      if (debug .eqv. .true.) then
        print *, "(", dSpaceOrbI1, dSpaceOrbIp2, "   | ",dSpaceOrbI2, dSpaceOrbIp1, ") = ", KVee(im_iq,in_ip)
        !tmpVal2 = tmpVal2 - KVee(im_iq, in_ip)
      endif
    endif

  end subroutine matrixElementsTwo

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

  subroutine readDipole(filePath, KDx, KDy, KDz)
    character(len=*), intent(in)        :: filePath
    real,             intent(out)     :: KDx(:,:), KDy(:,:), KDz(:,:)
    integer ::  i,j, xfile, yfile, zfile
   
    xfile = 100
    yfile = 101
    zfile = 102
    open(xfile,file=trim(filePath)//'x', status='old')
    open(yfile,file=trim(filePath)//'y', status='old')
    open(zfile,file=trim(filePath)//'z', status='old')
    do i=1,3
      read(xfile,*)
      read(yfile,*)
      read(zfile,*)
    enddo
    do i=1,NeleSpatialOrb
      read(xfile,*) (KDx(i,j),j=1,NeleSpatialOrb)
      read(yfile,*) (KDy(i,j),j=1,NeleSpatialOrb)
      read(zfile,*) (KDz(i,j),j=1,NeleSpatialOrb)
    enddo
    close(xfile)
    close(yfile)
    close(zfile)

    if (debug) then
      print *, "KDx"
      do i=1,NeleSpatialOrb
        print *, KDx(i,:)
      enddo
      print *, "KDy"
      do i=1,NeleSpatialOrb
        print *, KDy(i,:)
      enddo
      print *, "KDz"
      do i=1,NeleSpatialOrb
        print *, KDz(i,:)
      enddo
    endif


  end subroutine readDipole

end module operators
