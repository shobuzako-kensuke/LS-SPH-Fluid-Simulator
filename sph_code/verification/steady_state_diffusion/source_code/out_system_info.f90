subroutine out_system_info
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
    character(len=999) :: path_name, file_name

    !=========================!
    !  mkdir                  !
    !=========================!
    path_name = '../output/'//trim(adjustl(save_name))//'/system_info/'
    call mkdir(trim(adjustl(path_name)))

    !=========================!
    !  write                  !
    !=========================!
    file_name = trim(adjustl(path_name))//'system_info_NAMES.dat'
    open(10, file=file_name, status='replace', form='formatted')  ! save variable names as ascii

    file_name = trim(adjustl(path_name))//'system_info_VALUES.dat'
    open(11, file=file_name, status='replace', form='formatted')  ! save variable values as ascii

    ! write variable names
    write(10, *) 'write_step'
    write(10, *) 'Nx'
    write(10, *) 'Nh'
    write(10, *) 'rho_ref'
    write(10, *) 'x_rand'
    write(10, *) 'threshold'
    
    write(10, *) 'D'
    write(10, *) 'm'
    write(10, *) 'h'
    write(10, *) 'h_eff'
    write(10, *) 'WL_thick'
    write(10, *) 'dt'
    write(10, *) 'N_WL'
    write(10, *) 'N_inn'
    write(10, *) 'N_out'
    write(10, *) 'N_sys'

    ! write variable values in the same order as variable names
    write(11, *) write_step
    write(11, *) Nx
    write(11, *) Nh
    write(11, *) rho_ref
    write(11, *) x_rand
    write(11, *) threshold
    
    write(11, *) D
    write(11, *) m
    write(11, *) h
    write(11, *) h_eff
    write(11, *) WL_thick
    write(11, *) dt
    write(11, *) N_WL
    write(11, *) N_inn
    write(11, *) N_out
    write(11, *) N_sys

    close(10)
    close(11)

end subroutine out_system_info

! END !
