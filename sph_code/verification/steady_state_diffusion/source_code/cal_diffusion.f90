subroutine cal_diffusion
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
    integer :: RK_count, me, you, my_cell, your_cell, link_cell, num, N_link
    real(8) :: pi=acos(-1.0d0), x(N_you,3), W(N_you), r(N_you), f(N_you), Lap
    
    SP_f_now(:) = SP_f(:N_inn)  ! save function value at n step
    RK = 0.0d0  ! zero clear

    do RK_count = 1, 2  ! using 2nd-order Runge-Kutta method
        !=========================!
        !  cal diffusion          !
        !=========================!
        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(N_threads)
        !$OMP parallel default(none) &
        !$OMP private(me, you, my_cell, your_cell, link_cell, num, N_link) &
        !$OMP private(x, W, r, f, Lap) &
        !$OMP shared(N_inn, cell, cell_max, cell_9, N_cell, cell_info, m, h, N_you) &
        !$OMP shared(SP_x, SP_r, SP_f, SP_f_next, SP_f_now, pi, dt, RK, RK_count, WL_thick)
        !$OMP do
        do me = 1, N_inn  ! loop for inner particles

            !=========================!
            !  zero clear             !
            !=========================!
            N_link = 0
            x = 0.0d0
            W = 0.0d0
            r = 0.0d0
            f = 0.0d0

            my_cell = cell(me)  ! my cell

            do link_cell = 1, 9
                your_cell = cell_9(link_cell, my_cell)  ! neighboring 9 cells

                if ((your_cell < 1) .or. (cell_max < your_cell)) cycle  ! out of range
                if (N_cell(your_cell) == 0) cycle                       ! no particle in your_cell

                do num = 1, N_cell(your_cell)
                    you = cell_info(num, your_cell)  ! your name

                    if (me == you) cycle  ! if me=you, skip

                    N_link = N_link + 1   ! number of linked neighbors

                    if (N_link > N_you) stop ' [error] N_link > N_you @ cal_diffusion.f90'

                    x(N_link,1:2) = SP_x(me,:) - SP_x(you,:)  ! x_ij = x_i - x_j
                    x(N_link,3  ) = sqrt(x(N_link,1)**2.0d0 + x(N_link,2)**2.0d0)
#if defined(classical_SPH_Laplacian)
                    call cal_dW(x(N_link,3), h, W(N_link))    ! from lib_kernel_functions
#elif defined(LS_SPH_2ND) || defined(LS_SPH_3RD)
                    call cal_W(x(N_link,3), h, W(N_link))     ! from lib_kernel_functions
#endif
                    r(N_link) = SP_r(you)                     ! density (rho)
                    f(N_link) = SP_f(you) - SP_f(me)          ! function value
                enddo
            enddo


            !=========================!
            !  SPH approximation      !
            !=========================!
            if (N_link == 0) cycle  ! if there is no neighbors, skip

#if defined(classical_SPH_Laplacian)
            call cal_classical_Laplacian(N_you, N_link, m, x, W, r, f, Lap)
#elif defined(LS_SPH_2ND)
            call cal_LSSPH_typeA(5, N_you, N_link, m, h, x, W, r, f, Lap)
#elif defined(LS_SPH_3RD)
            call cal_LSSPH_typeA(9, N_you, N_link, m, h, x, W, r, f, Lap)
#endif

            !=========================!
            !  cal right-hand side    !
            !=========================!
            RK(me,RK_count) = dt* (Lap + 2.0d0*(pi**2.0d0)* sin(pi*(SP_x(me,1)-WL_thick)) &
                                                          * cos(pi*(SP_x(me,2)-WL_thick)))
        enddo
        !$OMP enddo
        !$OMP end parallel

        !=========================!
        !  cal next step          !
        !=========================!
        if (RK_count==1) then
            SP_f(:N_inn) = SP_f(:N_inn) + RK(:,1)
            call cal_VM_to_WL  ! renew wall
        elseif (RK_count==2) then
            SP_f_next(:) = SP_f_now(:) + 0.5d0* (RK(:,1) + RK(:,2))
            SP_f(:N_inn) = SP_f_now(:)
        endif
    enddo

end subroutine cal_diffusion

! END !