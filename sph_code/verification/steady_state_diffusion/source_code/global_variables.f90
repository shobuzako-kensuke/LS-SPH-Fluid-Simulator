module global_variables
    implicit none
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                             global variables                             !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    ! set_system.f90
    real(8), save :: D, m, h, h_eff, WL_thick, dt
    integer, save :: N_WL, N_inn, N_out, N_sys

    ! set_background_cell.f90
    integer, save :: cell_x_max, cell_y_max, cell_max, N_you

    ! main.f90
    integer, save :: step
    real(8), save :: start_time, tmp_time, residual_max


    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                              global arrays                               !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    ! set_system.f90
    integer, allocatable, save :: SP_kind(:)
    real(8), allocatable, save :: SP_x(:,:), SP_r(:), SP_f(:)
    real(8), allocatable, save :: SP_f_now(:), SP_f_next(:), RK(:,:), dif_f(:)

    ! set_VM.f90
    integer, allocatable, save :: VM_kind(:), VM_WL_name(:)
    real(8), allocatable, save :: VM_x(:,:), VM_f(:)
    
    ! set_background_cell.f90
    integer, allocatable, save :: cell(:), cell_VM(:), cell_9(:,:)
    integer, allocatable, save :: N_cell(:), N_cell_WL(:), cell_info(:,:), cell_info_fix(:,:)

    contains
    !=========================!
    !  deallocate             !
    !=========================!
    subroutine deallocate_array
        deallocate(SP_kind, SP_x, SP_r, SP_f, SP_f_now, SP_f_next, RK, dif_f)
        deallocate(VM_kind, VM_WL_name, VM_x, VM_f)
        deallocate(cell, cell_VM, cell_9, N_cell, N_cell_WL, cell_info, cell_info_fix)
    end subroutine deallocate_array

end module global_variables

! END !