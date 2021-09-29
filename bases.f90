module bases
  use params

  implicit none (type, external)
  external    ::  ssyevr
  private
  public      ::  diagSymMatrix

  contains
  subroutine diagSymMatrix(A, nEig, eigenVals, eigenVectors)
    real,                 intent(inout) ::  A(:,:)
    real,                 intent(inout) ::  eigenVectors(:,:)
    real,                 intent(inout) ::  eigenVals(:)
    integer,              intent(in)    ::  nEig
    integer                             ::  nEigFound
    real, dimension(:), allocatable     ::  work, work2
    integer, dimension(:), allocatable  ::  iwork, iwork2
    integer                             ::  dimA, lwork, liwork, &
                                            vl, vu, il, iu, info
    integer, dimension(:), allocatable  ::  isuppz
    real                                ::  abstol=-1.0

    il=1
    iu=nEig
    lwork=-1
    liwork=-1
    dimA = size(A,1)
    allocate(isuppz(2*dimA))
    allocate(work(1))
    allocate(iwork(1))

    !subroutine ssyevr	(	character 	JOBZ,
    !                     character 	RANGE,
    !                     character 	UPLO,
    !                     integer 	N,
    !                     real, dimension( lda, * ) 	A,
    !                     integer 	LDA,
    !                     real 	VL,
    !                     real 	VU,
    !                     integer 	IL,
    !                     integer 	IU,
    !                     real 	ABSTOL,
    !                     integer 	M,
    !                     real, dimension( * ) 	W,
    !                     real, dimension( ldz, * ) 	Z,
    !                     integer 	LDZ,
    !                     integer, dimension( * ) 	ISUPPZ,
    !                     real, dimension( * ) 	WORK,
    !                     integer 	LWORK,
    !                     integer, dimension( * ) 	IWORK,
    !                     integer 	LIWORK,
    !                     integer 	INFO 
    !                     )	

    call        ssyevr(   'V',            &
                          'I',            &
                          'U',            & 
                          dimA,           &
                          A,              &
                          dimA,           &
                          vl,             &
                          vu,             &
                          il,             &
                          iu,             &
                          abstol,         &
                          nEigFound,      &
                          eigenVals,      &
                          eigenVectors,   &
                          dimA,           &
                          isuppz,         &
                          work,           &
                          lwork,          &
                          iwork,          &
                          liwork,         &
                          info               )
    print *, "info: ", info
    lwork = work(1)
    deallocate(work)
    allocate(work2(lwork))
    liwork = iwork(1)
    deallocate(iwork)
    allocate(iwork2(liwork))
    
    call        ssyevr(   'V',            &
                          'I',            &
                          'U',            & 
                          dimA,           &
                          A,              &
                          dimA,           &
                          vl,             &
                          vu,             &
                          il,             &
                          iu,             &
                          abstol,         &
                          nEigFound,      &
                          eigenVals,      &
                          eigenVectors,   &
                          dimA,           &
                          isuppz,         &
                          work2,           &
                          lwork,          &
                          iwork2,          &
                          liwork,         &
                          info               )

    deallocate(iwork2)
    deallocate(work2)
    deallocate(isuppz)

    end subroutine diagSymMatrix

end module bases
