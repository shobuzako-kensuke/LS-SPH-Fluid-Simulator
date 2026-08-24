subroutine out_data(recording_ID)
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
    character(len=*), intent(in) :: recording_ID
    character(len=999) :: save_path, file_name

    !=========================!
    !  mkdir                  ! 
    !=========================!
    save_path = '../output/'//trim(adjustl(save_name))//'/data/'&
                //trim(adjustl(recording_ID))//'/'
    call mkdir(trim(adjustl(save_path)))

    !=========================!
    !  output                 ! 
    !=========================!
    if (trim(adjustl(recording_ID)) == '0') then              ! if the first loop
        file_name = trim(adjustl(save_path))//'VM_kind.dat'
        call write_binary(file_name, dble(VM_kind), N_VM,  1) ! VM_kind
        file_name = trim(adjustl(save_path))//'VM_xyz.dat'
        call write_binary(file_name, VM_xyz, N_VM, 2)         ! VM_xyz
    endif

    file_name = trim(adjustl(save_path))//'SP_kind.dat'
    call write_binary(file_name, dble(SP_kind), N_sys, 1)     ! SP_kind
    file_name = trim(adjustl(save_path))//'SP_xyz.dat'
    call write_binary(file_name, SP_xyz, N_sys, 2)            ! SP_xyz
    file_name = trim(adjustl(save_path))//'SP_uvw.dat'
    call write_binary(file_name, SP_uvw, N_sys, 2)            ! SP_uvw
    file_name = trim(adjustl(save_path))//'SP_rho.dat'
    call write_binary(file_name, SP_rho, N_sys, 1)            ! SP_rho
    file_name = trim(adjustl(save_path))//'SP_pre.dat'
    call write_binary(file_name, SP_pre, N_sys, 1)            ! SP_pre

#if defined(BOUSSINESQ_CONV)
    file_name = trim(adjustl(save_path))//'SP_tem.dat'
    call write_binary(file_name, SP_tem, N_sys, 1)            ! SP_tem
#endif

end subroutine out_data

! END !
