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
    write(10, *) 'start_step'
    write(10, *) 'end_step'
    write(10, *) 'write_step'
    write(10, *) 'total_step'

    write(10, *) 'N_sys'
    write(10, *) 'N_inn'
    write(10, *) 'N_out'
    write(10, *) 'N_VM'
    
    write(10, *) 'SP_mass'
    write(10, *) 'Delta_x'
    write(10, *) 'h'
    write(10, *) 'h_eff'
    write(10, *) 'WL_thick'
    write(10, *) 'W_ave'
    ! 14 items above

    write(10, *) 'zeta'
    write(10, *) 'xi'
    write(10, *) 'N_x'
    write(10, *) 'N_y'
    write(10, *) 'N_h'
    write(10, *) 'width'
    write(10, *) 'height'

    write(10, *) 'Delta_t_CFL'
    write(10, *) 'Delta_t_vis'
    write(10, *) 'Delta_t_th'
    write(10, *) 'Delta_t_CFL_relax'
    write(10, *) 'Delta_t_vis_relax'
    write(10, *) 'Delta_t_eff'

    write(10, *) 'rho_ref'
    write(10, *) 'vis_ref'
    write(10, *) 'K_ref'
    write(10, *) 'k_th'
    write(10, *) 'c_p'
    write(10, *) 'alpha'

    write(10, *) 'U_top'
    write(10, *) 'T_top'
    write(10, *) 'T_bot'
    write(10, *) 'gravity'

    write(10, *) 'Delta_temp'
    write(10, *) 'temp_ave'    
    write(10, *) 'sound_speed'
    write(10, *) 'kinematic_vis'
    write(10, *) 'thermal_dif'

    write(10, *) 'Re'
    write(10, *) 'Ra'
    write(10, *) 'Pr'

    ! write variable values [note: the same order as variable names]
    write(11, *) start_step
    write(11, *) end_step
    write(11, *) write_step
    write(11, *) total_step

    write(11, *) N_sys
    write(11, *) N_inn
    write(11, *) N_out
    write(11, *) N_VM
    
    write(11, *) SP_mass
    write(11, *) Delta_x
    write(11, *) h
    write(11, *) h_eff
    write(11, *) WL_thick
    write(11, *) W_ave
    ! 14 items above

    write(11, *) zeta
    write(11, *) xi
    write(11, *) N_x
    write(11, *) N_y
    write(11, *) N_h
    write(11, *) width
    write(11, *) height

    write(11, *) Delta_t_CFL
    write(11, *) Delta_t_vis
    write(11, *) Delta_t_th
    write(11, *) Delta_t_CFL_relax
    write(11, *) Delta_t_vis_relax
    write(11, *) Delta_t_eff

    write(11, *) rho_ref
    write(11, *) vis_ref
    write(11, *) K_ref
    write(11, *) k_th
    write(11, *) c_p
    write(11, *) alpha

    write(11, *) U_top
    write(11, *) T_top
    write(11, *) T_bot
    write(11, *) gravity

    write(11, *) Delta_temp
    write(11, *) temp_ave
    write(11, *) sound_speed
    write(11, *) kinematic_vis
    write(11, *) thermal_dif

    write(11, *) Re
    write(11, *) Ra
    write(11, *) Pr

    close(10)
    close(11)

end subroutine out_system_info

! END !
