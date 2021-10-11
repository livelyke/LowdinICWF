module saving
  implicit none (type, external)
  external    ::  ssymv 
  private

  public    ::  output

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

end module saving
