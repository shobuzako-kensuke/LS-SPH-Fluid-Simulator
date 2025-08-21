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
    if (trim(adjustl(recording_ID)) == '0') then               ! if the first loop
        file_name = trim(adjustl(save_path))//'VM_kind.dat'
        call write_binary(file_name, dble(VM_kind), N_out, 1)  ! VM_kind
        file_name = trim(adjustl(save_path))//'VM_x.dat'
        call write_binary(file_name, VM_x, N_out, 2)           ! VM_xyz

        file_name = trim(adjustl(save_path))//'SP_kind.dat'
        call write_binary(file_name, dble(SP_kind), N_sys, 1)  ! SP_kind
        file_name = trim(adjustl(save_path))//'SP_x.dat'
        call write_binary(file_name, SP_x, N_sys, 2)           ! SP_xyz
        file_name = trim(adjustl(save_path))//'SP_r.dat'
        call write_binary(file_name, SP_r, N_sys, 1)           ! SP_r
    endif

    file_name = trim(adjustl(save_path))//'SP_f.dat'
    call write_binary(file_name, SP_f, N_sys, 1)               ! SP_f
end subroutine out_data

! END !
