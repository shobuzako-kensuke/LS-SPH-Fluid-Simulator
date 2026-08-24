subroutine set_background_cell
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
    integer :: me, cell_x, cell_y, my_cell, cell_cap
    real(8) :: x_max, y_max

    !=========================!
    !  cell max               !
    !=========================!
    x_max = maxval(SP_x(:,1))
    y_max = maxval(SP_x(:,2))
    cell_x_max = ceiling(x_max / h_eff)
    cell_y_max = ceiling(y_max / h_eff)
    cell_max   = cell_x_max + (cell_y_max - 1) * cell_x_max  ! max cell number

    !=========================!
    !  allocate               !
    !=========================!
    allocate(cell     (N_sys))       ! cell   (me) = my cell number
    allocate(cell_VM  (N_out))       ! cell_VM(me) = cell number of VM
    allocate(cell_9   (9,cell_max))  ! cell_9(9, my cell) = neighbor 9 cells
    allocate(N_cell   (cell_max))    ! N_cell   (my cell) = number of residents
    allocate(N_cell_WL(cell_max))    ! N_cell_WL(my_cell) = number of residents for WL

    !=========================!
    !  zero clear             !
    !=========================!
    cell      = 0
    cell_VM   = 0
    cell_9    = 0
    N_cell    = 0
    N_cell_WL = 0

    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me, cell_x, cell_y, my_cell) &
    !$OMP shared(N_sys, N_inn, N_out, SP_x, VM_x, h_eff, cell_x_max, cell_max, cell_cap, N_you) &
    !$OMP shared(N_cell, N_cell_WL, cell, cell_VM, cell_info, cell_info_fix, cell_9)
    !=========================!
    !  count for capacity     !
    !=========================!
    !$OMP do
    do me = 1, N_sys
        cell_x = ceiling(SP_x(me,1) / h_eff)
        cell_y = ceiling(SP_x(me,2) / h_eff)
        my_cell = cell_x + (cell_y - 1) * cell_x_max  ! my cell number
        cell(me) = my_cell                            ! store my cell number
    enddo
    !$OMP enddo
    !$OMP barrier

    !=========================!
    !  allocate cell_info     !
    !=========================!
    !$OMP single
    do me = 1, N_sys
        my_cell = cell(me)
        N_cell(my_cell) = N_cell(my_cell) + 1  ! count
    enddo

    cell_cap = maxval(N_cell(:))  ! capacity in a cell
    N_you = 12*cell_cap           ! maximum number of neighbors in a cell

    allocate(cell_info    (2*cell_cap, cell_max))  ! cell_info    (names, my_cell)
    allocate(cell_info_fix(2*cell_cap, cell_max))  ! cell_info_fix(names, my_cell)    
    cell_info     = 0  ! zero clear
    cell_info_fix = 0
    
    !=========================!
    !  cell_info_fix          !
    !=========================!
    do me = N_inn+1, N_sys  ! loop for wall (outer) particles
        my_cell = cell(me)
        N_cell_WL(my_cell) = N_cell_WL(my_cell) + 1  ! count
        cell_info_fix(N_cell_WL(my_cell), my_cell) = me
    enddo
    !$OMP end single
    !$OMP barrier

    !=========================!
    !  cell_9                 !
    !=========================!
    !$OMP do
    do my_cell = 1, cell_max
        cell_9(1, my_cell) = my_cell - cell_x_max - 1 ! lower left
        cell_9(2, my_cell) = my_cell - cell_x_max     ! lower center
        cell_9(3, my_cell) = my_cell - cell_x_max + 1 ! lower right
        cell_9(4, my_cell) = my_cell - 1              ! left
        cell_9(5, my_cell) = my_cell                  ! center
        cell_9(6, my_cell) = my_cell + 1              ! right
        cell_9(7, my_cell) = my_cell + cell_x_max - 1 ! upper left
        cell_9(8, my_cell) = my_cell + cell_x_max     ! upper center
        cell_9(9, my_cell) = my_cell + cell_x_max + 1 ! upper right
    enddo
    !$OMP enddo
    !$OMP barrier

    !=========================!
    !  cell_VM                !
    !=========================!
    !$OMP do
    do me = 1, N_out
        cell_x = ceiling(VM_x(me,1) / h_eff)
        cell_y = ceiling(VM_x(me,2) / h_eff)
        my_cell = cell_x + (cell_y - 1) * cell_x_max  ! my cell number
        cell_VM(me) = my_cell                         ! store my cell number
    enddo
    !$OMP enddo
    !$OMP end parallel
  
end subroutine set_background_cell

! END !