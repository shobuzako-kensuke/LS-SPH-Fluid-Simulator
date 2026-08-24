module global_variables
    implicit none
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                             global variables                             ! 
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !  set_system_info.f90
    integer, save :: N_WL, total_step                        ! SPH info
    real(8), save :: Delta_x, SP_mass, h, h_eff, WL_thick    ! SPH info
    real(8), save :: Delta_temp, temp_ave, sound_speed       ! physical info
    real(8), save :: kinematic_vis, thermal_dif              ! physical info
    real(8), save :: Re, Ra, Pr                              ! physical info
    real(8), save :: Delta_t_CFL, Delta_t_vis, Delta_t_th    ! time step
    real(8), save :: Delta_t_CFL_relax, Delta_t_vis_relax    ! time step
    real(8), save :: Delta_t_eff                             ! time step

    !  set_xyz_kind.f90
    integer, save :: N_sys, N_inn, N_out
    
    !  set_VM_xy_kind.f90
    integer, save :: N_VM

    !  set_physical_quantity.f90
    real(8), save :: W_ave

    !  set_background_cell.f90
    integer, save :: cell_x_max, cell_y_max, cell_max, cell_cap, N_you

    !  main.f90
    integer, save :: step
    real(8), save :: start_time, tmp_time


    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                              global arrays                               ! 
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !  set_xy_kind.f90
    integer, allocatable, save :: SP_kind(:)
    real(8), allocatable, save :: SP_xyz(:,:)

    !  set_VM_xy_kind.f90
    integer, allocatable, save :: VM_kind(:), VM_WL_name(:)
    real(8), allocatable, save :: VM_xyz(:,:), VM_rho(:), VM_pre(:), VM_uvw(:,:), VM_tem(:)

    !  set_physical_quantity.f90
    real(8), allocatable, save :: SP_rho(:), SP_pre(:), SP_uvw(:,:), SP_tem(:)
    
    !  set_background_cell.f90
    integer, allocatable, save :: cell(:), cell_VM(:), cell_9(:,:)
    integer, allocatable, save :: cell_info(:,:), cell_info_fix(:,:)
    integer, allocatable, save :: N_cell(:), N_cell_WL(:)


    contains
    !=========================!
    !  subroutine deallocate  ! 
    !=========================!
    subroutine deallocate_array
        deallocate(SP_kind)
        deallocate(SP_xyz)
        deallocate(VM_kind, VM_WL_name)
        deallocate(VM_xyz, VM_rho, VM_pre, VM_uvw, VM_tem)
        deallocate(SP_rho, SP_pre, SP_uvw, SP_tem)
        deallocate(cell, cell_VM, cell_9)
        deallocate(cell_info, cell_info_fix)
        deallocate(N_cell, N_cell_WL)
    end subroutine deallocate_array

end module global_variables

! END !