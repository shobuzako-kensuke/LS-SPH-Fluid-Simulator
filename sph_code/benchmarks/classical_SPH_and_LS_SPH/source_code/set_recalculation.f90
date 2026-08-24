subroutine set_recalculation
    !=========================!
    !  module                 ! 
    !=========================!
    use input
    use global_variables
    use lib_file_operations

    !=========================!
    !  local variables        ! 
    !=========================!
    implicit none
    character(len=999) :: read_path, file_name, file_ID
    integer :: i
    real(8) :: tmp
    real(8), allocatable :: SP_kind_tmp(:)

    !=========================!
    !  read settings          !
    !=========================!
    file_name = '../output/'//trim(adjustl(save_name))//'/system_info/system_info_VALUES.dat'
    open(10, file=file_name, form='formatted', status='old')

    do i = 1, 14
        read(10,*) tmp ! read system_info.dat

        if (i == 5) then
            N_sys = int(tmp)
        endif

        if (i == 6) then
            N_inn = int(tmp)
        endif

        if (i == 7) then
            N_out = int(tmp)
        endif

        if (i == 14) then
            W_ave = tmp
        endif
    enddo
    close(10)

    !=========================!
    !  allocate               !
    !=========================!
    allocate(SP_xyz     (N_sys, 2))
    allocate(SP_kind    (N_sys))
    allocate(SP_kind_tmp(N_sys))  ! for read real array
    allocate(SP_rho     (N_sys))
    allocate(SP_pre     (N_sys))
    allocate(SP_uvw     (N_sys,2))
    allocate(SP_tem     (N_sys))

    SP_xyz  = 0.0d0 ! zero clear
    SP_kind = 0
    SP_rho  = 0.0d0
    SP_pre  = 0.0d0
    SP_uvw  = 0.0d0
    SP_tem  = 0.0d0

    !=========================!
    !  read files             !
    !=========================!
    ! read last data
    write(file_ID,*) start_step
    read_path = '../output/'//trim(adjustl(save_name))//'/data/'//trim(adjustl(file_ID))//'/'

    ! SP_kind
    file_name = trim(adjustl(read_path))//'SP_kind.dat'
    call read_binary(file_name, SP_kind_tmp, N_sys, 1)
    SP_kind(:) = int(SP_kind_tmp(:))  ! float >> integer
    deallocate(SP_kind_tmp)

    ! SP_xyz
    file_name = trim(adjustl(read_path))//'SP_xyz.dat'
    call read_binary(file_name, SP_xyz, N_sys, 2)

    ! SP_uvw
    file_name = trim(adjustl(read_path))//'SP_uvw.dat'
    call read_binary(file_name, SP_uvw, N_sys, 2)

    ! SP_rho
    file_name = trim(adjustl(read_path))//'SP_rho.dat'
    call read_binary(file_name, SP_rho, N_sys, 1)

    ! SP_pre
    file_name = trim(adjustl(read_path))//'SP_pre.dat'
    call read_binary(file_name, SP_pre, N_sys, 1)

    ! SP_tem
#if defined(TEMP_ON)
    file_name = trim(adjustl(read_path))//'SP_tem.dat'
    call read_binary(file_name, SP_tem, N_sys, 1)
#endif

end subroutine set_recalculation

! END !