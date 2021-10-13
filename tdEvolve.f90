module tdEvolve
  use params,     only : NSlater,NnucOrb, dtImag, dt, electronicOnly, &
                         nAxis, debug, finalTime, kick, kickDir, kappa
  use bases,      only : diagSymMatrix
  use saving,     only : output
  use operators,  only : MDx, MDy, MDz
  implicit none (type, external)
  external    ::  ssymv 
  private
  public  ::  tdEvolveImag, tdEvolveReal

  real, public, allocatable, dimension(:,:) :: C

  contains 

  subroutine tdEvolveReal(H,C)
    real, intent(inout)   :: C(:,:)
    real, intent(in)      :: H(:,:)
    real, dimension(:,:), allocatable :: CtmpR, CtmpI
    real        :: time
    integer     :: n, ti, saveInterval

    n = size(H,1)
    allocate(CtmpR(n,4))
    CtmpR = 0
    allocate(CtmpI(n,4))
    CtmpI = 0

    time=0
    ti=0
    saveInterval = int(0.5/dt)
    call output(time, H, C, MDx, MDy, MDz, CtmpR, CtmpI)
    if ( kick .eqv. .true.) then
      if (kickDir == 1) then
        call ssymv('U',n, -1.0*kappa, MDx,n,C(:,1) ,1,0.0,CtmpI(:,1),1)
      elseif (kickDir == 2) then
        call ssymv('U',n, -1.0*kappa, MDy,n,C(:,1) ,1,0.0,CtmpI(:,1),1)
      elseif (kickDir == 3) then
        call ssymv('U',n, -1.0*kappa, MDz,n,C(:,1) ,1,0.0,CtmpI(:,1),1)
      endif
      C(:,2) = CtmpI(:,1)
    endif
        
    do while (time <= finalTime + dt)
      call rk4(H, C, dt, CtmpR, CtmpI)
      time = time + dt
      ti = ti+1
      if ( mod(ti,saveInterval) == 0 ) then
        call output(time, H, C, MDx, MDy, MDz, CtmpR, CtmpI)
      endif
    enddo
    

  end subroutine tdEvolveReal


  subroutine tdEvolveImag(H)
    real, intent(in)  ::  H(:,:)
    integer           ::  n, iter, nEig
    real, dimension(:),allocatable  :: HC, eigenVals
    real, dimension(:,:),allocatable  :: eigenVectors, Htmp
    real              ::  Enew,Eold

    print *, "Diagonalizing Hamiltonian"
    n=size(H,1)
    allocate(C(n,2))
    C = 0
    allocate(HC(n))
    HC=0
    allocate(Htmp(n,n))
    Htmp=0
    Htmp = H
    !call random_number(C)
    !C = (/ (rand(2*i), i=1,n) /)
    !C = C/sqrt(sum(C*C))

    !call ssymv('u',n,1.0,H,n,C,1,0.0,HC,1)
    !Eold = sum(C*HC)
    !Enew = Eold + 1
    !iter = 0
    !do while (abs(Eold-Enew) > 1e-12 .and. iter<=1e6)
    !  Eold = Enew 
    !  call rk4Imag(H, C,dtImag)
    !  call ssymv('u',n,1.0,H,n,C,1,0.0,HC,1)
    !  Enew = sum(C*HC)
    !  iter = iter+1
    !  if(mod(iter,100)==0) then
    !   print *, "iter: ", iter, "E:", Enew
    !  endif
    !enddo
    !print *, "iter: ", iter, "E:", Enew

    nEig = min(20,n)
    allocate(eigenVals(nEig))
    eigenVals = 0
    allocate(eigenVectors(n,nEig))
    eigenVectors = 0
  
    if (debug .eqv. .true.) then
      print *, "H"
      do iter=1,n
        print *, Htmp(iter,:)
      enddo
    endif
    call diagSymMatrix(Htmp,nEig,eigenVals,eigenVectors)
    if (electronicOnly .eqv. .true.) then
      print *, nAxis(1), eigenVals
    else
      print *, eigenVals
    endif

    C(:,1) = eigenVectors(:,1)/sum(eigenVectors(:,1)*eigenVectors(:,1))
    deallocate(eigenVals)
    deallocate(eigenVectors)
 
  end subroutine

  subroutine rk4Imag(H, C, Ctmp, dt)
    real, intent(inout)   :: C(:,:), Ctmp(:,:)
    real, intent(in)      :: dt, H(:,:)
    integer               :: i, n
    n = size(C,1)
    Ctmp = 0

    !> y := Ax + \beta*y
    !> ssymv ( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
    !do i=1,4
    call ssymv('U',n,1.0,H,n,C(:,1)                     ,1,0.0,Ctmp(:,1),1)
    call ssymv('U',n,1.0,H,n,C(:,1) + (dt/2.0)*Ctmp(:,1),1,0.0,Ctmp(:,2),1)
    call ssymv('U',n,1.0,H,n,C(:,1) + (dt/2.0)*Ctmp(:,2),1,0.0,Ctmp(:,3),1)
    call ssymv('U',n,1.0,H,n,C(:,1) + (dt)    *Ctmp(:,2),1,0.0,Ctmp(:,4),1)

    C(:,1) = C(:,1) - (dt/6.0)*(Ctmp(:,1) + 2*Ctmp(:,2) + 2*Ctmp(:,3) + Ctmp(:,4))

    C(:,1) = C(:,1)/sqrt(sum(C(:,1)*C(:,1)))

  end subroutine
  
  subroutine rk4(H, C, dt, CtmpR, CtmpI)
    real, intent(inout)   :: C(:,:)
    real, intent(inout)   :: CtmpR(:,:), CtmpI(:,:)
    real, intent(in)      :: H(:,:)
    real, intent(in)      :: dt
    integer               :: i, n

    n = size(C,1)

    CtmpR = 0
    CtmpI = 0

    !> Split up the real and imaginary parts
    !> \dot{C} = -iHC = -iHC_r + HC_i
    !> y := Ax + \beta*y
    !> ssymv ( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
    !do i=1,4
    call ssymv('U',n,-1.0,H,n,C(:,1)                      ,1,0.0,CtmpI(:,1),1)
    call ssymv('U',n, 1.0,H,n,C(:,2)                      ,1,0.0,CtmpR(:,1),1)

    call ssymv('U',n,-1.0,H,n,C(:,1) + (dt/2.0)*CtmpR(:,1),1,0.0,CtmpI(:,2),1)
    call ssymv('U',n, 1.0,H,n,C(:,2) + (dt/2.0)*CtmpI(:,1),1,0.0,CtmpR(:,2),1)

    call ssymv('U',n,-1.0,H,n,C(:,1) + (dt/2.0)*CtmpR(:,2),1,0.0,CtmpI(:,3),1)
    call ssymv('U',n, 1.0,H,n,C(:,2) + (dt/2.0)*CtmpI(:,2),1,0.0,CtmpR(:,3),1)

    call ssymv('U',n,-1.0,H,n,C(:,1) + (dt)   *CtmpR(:,3),1,0.0,CtmpI(:,4),1)
    call ssymv('U',n, 1.0,H,n,C(:,2) + (dt   )*CtmpI(:,3),1,0.0,CtmpR(:,4),1)

    C(:,1) = C(:,1) + (dt/6.0)*(CtmpR(:,1) + 2*CtmpR(:,2) + 2*CtmpR(:,3) + CtmpR(:,4))
    C(:,2) = C(:,2) + (dt/6.0)*(CtmpI(:,1) + 2*CtmpI(:,2) + 2*CtmpI(:,3) + CtmpI(:,4))

  end subroutine rk4
   

end module tdEvolve
