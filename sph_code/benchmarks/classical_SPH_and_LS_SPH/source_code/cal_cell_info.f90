subroutine cal_cell_info
    !=========================!
    !  module                 ! 
    !=========================!
    !$use omp_lib
    use input
    use global_variables

    !=========================!
    !  local variables        ! 
    !=========================!
    implicit none
    integer :: me, cell_x, cell_y, my_cell, start_point, end_point

    !=========================!
    !  zero clear             ! 
    !=========================!
    cell(1:N_inn) = 0 ! clear ONLY inner particles 
    N_cell = 0
    cell_info = 0

    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me, cell_x, cell_y, my_cell, start_point, end_point) &
    !$OMP shared(N_inn, h_eff, cell_x_max, cell_max, cell_info, cell_info_fix) &
    !$OMP shared(SP_xyz, cell, N_cell, N_cell_WL, SP_kind)

    !=========================!
    !  cal my cell            !
    !=========================!
    !$OMP do
    do me = 1, N_inn ! ONLY inner particles
        cell_x = ceiling(SP_xyz(me,1) / h_eff)
        cell_y = ceiling(SP_xyz(me,2) / h_eff)
        my_cell = cell_x + (cell_y - 1) * cell_x_max  ! my cell number
        cell(me) = my_cell                            ! store my cell number
    enddo
    !$OMP enddo
    !$OMP barrier

    !=========================!
    !  cell_info              !
    !=========================!
    !$OMP single
    do me = 1, N_inn
        my_cell = cell(me)
        N_cell(my_cell) = N_cell(my_cell) + 1      ! count
        cell_info(N_cell(my_cell), my_cell) = me   ! make cell_info
    enddo
    !$OMP end single
    !$OMP barrier

    !=========================!
    !  combine cell_info_fix  !
    !=========================!
    !$OMP do
    do my_cell = 1, cell_max
        if (N_cell_WL(my_cell) > 0) then  ! if my_cell includes wall particles
            start_point = N_cell(my_cell) + 1
            end_point   = N_cell(my_cell) + N_cell_WL(my_cell)
            cell_info(start_point:end_point, my_cell) = cell_info_fix(1:N_cell_WL(my_cell), my_cell)
            N_cell(my_cell) = end_point   ! sum N_cell
        endif
    enddo
    !$OMP enddo
    !$OMP end parallel

end subroutine cal_cell_info

! END !