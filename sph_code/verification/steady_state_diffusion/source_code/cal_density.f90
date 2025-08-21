subroutine cal_density
    !=========================!
    !  module                 !
    !=========================!
    !$use omp_lib
    use input
    use global_variables
    use lib_kernel_functions

    !=========================!
    !  local variables        !
    !=========================!
    implicit none
    integer :: me, you, my_cell, your_cell, link_cell, num, N_link
    real(8) :: x(N_you,3), W(N_you), r(N_you)
    real(8), allocatable :: cal_r(:)
    allocate(cal_r(N_inn))
    cal_r = 0.0d0      ! zero clear

    SP_r(:) = rho_ref  ! reference density (inner + outer particles)

#if defined(VARIABLE_DENSITY)
    !=========================!
    !  cal density            !
    !=========================!
    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me, you, my_cell, your_cell, link_cell, num, N_link) &
    !$OMP private(x, W, r) &
    !$OMP shared(N_inn, cell, cell_max, cell_9, N_cell, cell_info, m, h, N_you) &
    !$OMP shared(SP_x, SP_r, cal_r)
    !$OMP do
    do me = 1, N_inn  ! loop for inner particles

        !=========================!
        !  zero clear             !
        !=========================!
        N_link = 0
        x = 0.0d0
        W = 0.0d0
        r = 0.0d0

        my_cell = cell(me)  ! my cell

        do link_cell = 1, 9
            your_cell = cell_9(link_cell, my_cell)  ! neighboring 9 cells

            if ((your_cell < 1) .or. (cell_max < your_cell)) cycle  ! out of range
            if (N_cell(your_cell) == 0) cycle                       ! no particle in your_cell

            do num = 1, N_cell(your_cell)
                you = cell_info(num, your_cell)  ! your name

                ! if (me == you) cycle  ! if me=you, skip

                N_link = N_link + 1  ! number of linked neighbors

                if (N_link > N_you) stop ' [error] N_link > N_you @ cal_density.f90'

                x(N_link,1:2) = SP_x(me,:) - SP_x(you,:)  ! x_ij = x_i - x_j
                x(N_link,3  ) = sqrt(x(N_link,1)**2.0d0 + x(N_link,2)**2.0d0)
                call cal_W(x(N_link,3), h, W(N_link))     ! from lib_kernel_functions
                r(N_link) = SP_r(you)                     ! density (rho)
            enddo
        enddo

        !=========================!
        !  SPH approximation      !
        !=========================!
        if (N_link == 0) cycle  ! if there is no neighbors, skip
        
        call cal_CSPH_0th(N_you, N_link, m, W, r, r, cal_r(me))
    enddo
    !$OMP enddo
    !$OMP end parallel

    SP_r(:N_inn) = cal_r(:)  ! substitution only for inner particles
#endif

    deallocate(cal_r)
    
end subroutine cal_density

! END !