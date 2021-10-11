module tdEvolve
  use params,     only : NSlater,NnucOrb, dtImag, electronicOnly, nAxis, debug
  use operators,  only : H
  use bases,      only : diagSymMatrix
  implicit none (type, external)
  external    ::  ssymv 
  private
  public  ::  tdEvolveImag

  real, public, allocatable, dimension(:) :: C

  contains 

  subroutine tdEvolveImag(H)
    real, intent(in)  ::  H(:,:)
    integer           ::  n, iter, nEig
    real, dimension(:),allocatable  :: HC, eigenVals
    real, dimension(:,:),allocatable  :: eigenVectors, Htmp
    real              ::  Enew,Eold

    n=size(H,1)
    allocate(C(n))
    C=0
    allocate(HC(n))
    HC=0
    allocate(Htmp(n,n))
    Htmp=0
    Htmp = H
    call random_number(C)
    !C = (/ (rand(2*i), i=1,n) /)
    C = C/sqrt(sum(C*C))

    !call ssymv('u',n,1.0,H,n,C,1,0.0,HC,1)
    !Eold = sum(C*HC)
    !Enew = Eold + 1
    !iter = 0
    !do while (abs(Eold-Enew) > 1e-12 .and. iter<=1e6)
    !  Eold = Enew 
    !  call rk4Imag(C,dtImag)
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
    endif
    deallocate(eigenVals)
    deallocate(eigenVectors)
 
  end subroutine

  subroutine rk4Imag(C,dt)
    real, intent(inout)   :: C(:)
    real, intent(in)      :: dt
    real, dimension(:,:), allocatable :: Ctmp
    integer               :: i, n
    n = size(C)
    allocate(Ctmp(n,4))
    Ctmp = 0

    !> y := Ax + \beta*y
    !> ssymv ( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
    !do i=1,4
    call ssymv('U',n,1.0,H,n,C                     ,1,0.0,Ctmp(:,1),1)
    call ssymv('U',n,1.0,H,n,C + (dt/2.0)*Ctmp(:,1),1,0.0,Ctmp(:,2),1)
    call ssymv('U',n,1.0,H,n,C + (dt/2.0)*Ctmp(:,2),1,0.0,Ctmp(:,3),1)
    call ssymv('U',n,1.0,H,n,C + (dt)    *Ctmp(:,2),1,0.0,Ctmp(:,4),1)

    C = C - (dt/6.0)*(Ctmp(:,1) + 2*Ctmp(:,2) + 2*Ctmp(:,3) + Ctmp(:,4))

    C = C/sqrt(sum(C*C))
    deallocate(Ctmp)
  end subroutine
     

end module tdEvolve
