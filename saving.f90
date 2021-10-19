module saving
  implicit none (type, external)
  external    ::  ssymv 
  private

  public    ::  output, outputNucDen

  contains

  subroutine output(time, H, C, MDx, MDy, MDz, CtmpR, CtmpI)
    real, intent(in)    :: time
    real, intent(in)    :: H(:,:), C(:,:), MDx(:,:), MDy(:,:), MDz(:,:)
    real, intent(inout) :: CtmpR(:,:), CtmpI(:,:) !< Recycle these for tmp storage
    logical             :: dirExist
    integer             :: dipoleFile, energyFile, n, i
    real                :: tmp1, tmp2, tmp3
    
    if ( time==0 ) then
      inquire (file='output_icwf/dipole.dat', exist=dirExist)
      if (.not. dirExist) then
        call execute_command_line('mkdir -p output_icwf')
        call execute_command_line('touch output_icwf/dipole.dat')
      endif
      inquire (file='output_icwf/energy.dat', exist=dirExist)
      if (.not. dirExist) then
        call execute_command_line('touch output_icwf/energy.dat')
      endif
    endif

    dipoleFile = 100
    energyFile = 101

    open(dipoleFile,file='output_icwf/dipole.dat',status='old')
    open(energyFile,file='output_icwf/energy.dat',status='old')

    n=size(C,1)

    !> get C'*H*C
    !> (C_r - iC_i)*H*(C_r + iC_i)
    !> H*C_r
    call ssymv('U',n,1.0,H,n,C(:,1)                     ,1,0.0,CtmpR(:,1),1)
    !> C_r*(H*C_r)
    tmp1 =  sum(C(:,1)*CtmpR(:,1))
    !> H*C_i
    call ssymv('U',n,1.0,H,n,C(:,2)                     ,1,0.0,CtmpI(:,1),1)
    !> C_i*(H*C_i)
    tmp1 =  tmp1 + sum(C(:,2)*CtmpI(:,1))

    tmp2 = sum(C(:,1)**2) + sum(C(:,2)**2)

    write(energyFile, *) time, tmp1, tmp2
    
    !> Get Dipole Moments
    !> (C_r - iC_i)*MDx*(C_r + iC_i)
    call ssymv('U',n,1.0,MDx,n,C(:,1)     ,1,0.0,CtmpR(:,2),1)
    !> C_r*(MDx*C_r)
    tmp1 =  sum(C(:,1)*CtmpR(:,2))
    !> H*C_i
    call ssymv('U',n,1.0,MDx,n,C(:,2)     ,1,0.0,CtmpI(:,2),1)
    !> C_i*(MDx*C_i)
    tmp1 =  tmp1 + sum(C(:,2)*CtmpI(:,2))
    
    call ssymv('U',n,1.0,MDy,n,C(:,1)     ,1,0.0,CtmpR(:,3),1)
    !> C_r*(MDy*C_r)
    tmp2 =  sum(C(:,1)*CtmpR(:,3))
    !> H*C_i
    call ssymv('U',n,1.0,MDy,n,C(:,2)     ,1,0.0,CtmpI(:,3),1)
    !> C_i*(MDy*C_i)
    tmp2 =  tmp2 + sum(C(:,2)*CtmpI(:,3))

    call ssymv('U',n,1.0,MDz,n,C(:,1)     ,1,0.0,CtmpR(:,4),1)
    !> C_r*(MDz*C_r)
    tmp3 =  sum(C(:,1)*CtmpR(:,4))
    !> H*C_i
    call ssymv('U',n,1.0,MDz,n,C(:,2)     ,1,0.0,CtmpI(:,4),1)
    !> C_i*(MDz*C_i)
    tmp3 =  tmp3 + sum(C(:,2)*CtmpI(:,4))

    
    write(dipoleFile,*) time, tmp1, tmp2, tmp3 


  end subroutine output

  subroutine outputNucDen(time,C,chi, nAxis)
    real, intent(in)  ::  C(:,:), chi(:,:), time, nAxis(:)
    real, dimension(:), allocatable :: nucDen
    integer, dimension(:), allocatable :: nucRange, slaterIJ, slaterIJp
    integer               ::  i,j,jp, NSlater, Nnuc, nDim
    logical               :: dirExist

    Nnuc    = size(chi,2)
    NSlater = size(C,1)/Nnuc
    nDim    = size(chi,1)

    allocate(nucDen(nDim))
    nucDen = 0

    allocate(nucRange(Nnuc))
    nucRange = (/ (i, i=1,Nnuc) /)
    allocate(slaterIJ(Nnuc))
    slaterIJ = 0
    allocate(slaterIJp(Nnuc))
    slaterIJp = 0
    
    do i=1,NSlater
      slaterIJ = nucRange + Nnuc*(i-1)
      slaterIJp = nucRange + Nnuc*(i-1)
      
      !> C^*_{IJ'}C_{IJ}\chi_{J'}chi_J
      do j=1,Nnuc
        do jp=1,Nnuc
          nucDen = nucDen + sum(C(slaterIJp(jp),:)*C(slaterIJ(j),:)) &
                          *chi(:,j)*chi(:,jp)
        enddo
      enddo
    enddo

    if ( time==0 ) then
      inquire (file='output_icwf/nucDen.dat', exist=dirExist)
      if (.not. dirExist) then
        call execute_command_line('mkdir -p output_icwf')
        call execute_command_line('touch output_icwf/nucDen.dat')
      endif
      open(100,file='output_icwf/nucDen.dat', status='old')
      do i=1,Nnuc
        write(100,*) nAxis(i), nucDen(i)
      enddo
      close(100)
    else
      open(100,file='output_icwf/tmp.dat', status='new')
      do i=1,nDim
        write(100,*) nucDen(i)
      enddo
      call execute_command_line('paste output_icwf/nucDen.dat output_icwf/tmp.dat > output_icwf/nucDen.dat')
      close(100)
      call execute_command_line('rm output_icwf/tmp.dat')
      
    endif

  end subroutine outputNucDen

end module saving
