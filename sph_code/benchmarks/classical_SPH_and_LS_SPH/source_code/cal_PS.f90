subroutine cal_PS
    !=========================!
    !  module                 ! 
    !=========================!
    !$use omp_lib
    use input
    use global_variables
    use lib_kernel_function

    !=========================!
    !  local variables        ! 
    !=========================!
    implicit none
    integer :: me, you, my_cell, your_cell, link_cell, num, N_link
    real(8) :: x_link(N_you), y_link(N_you), r_link(N_you)
    real(8) :: W_link(N_you), dW_link(N_you)
    real(8) :: rho_link(N_you), uvw_link(2,N_you), tem_link(N_you)
    real(8) :: d_x(2), d_r(5), d_u(5), d_v(5), d_T(5)
    real(8) :: V_max, coe_PS
    real(8) :: V_tmp(N_inn), SP_xyz_new(N_inn,2), SP_uvw_new(N_inn,2)
    real(8) :: SP_rho_new(N_inn), SP_tem_new(N_inn)

    ! zero clear
    V_tmp = 0.0d0
    SP_xyz_new = 0.0d0
    SP_uvw_new = 0.0d0
    SP_rho_new = 0.0d0
    SP_tem_new = 0.0d0

    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me, you, my_cell, your_cell, link_cell, num, N_link) &
    !$OMP private(x_link, y_link, r_link, W_link, dW_link, rho_link, uvw_link, tem_link) &
    !$OMP private(d_x, d_r, d_u, d_v, d_T) &
    !$OMP shared(N_inn, cell, cell_max, cell_9, N_cell, cell_info, h, N_you) &
    !$OMP shared(SP_kind, SP_xyz, SP_rho, SP_pre, SP_uvw, SP_tem) &
    !$OMP shared(V_tmp, V_max, coe_PS, Delta_t_eff) &
    !$OMP shared(SP_xyz_new, SP_uvw_new, SP_rho_new, SP_tem_new)

    !=========================!
    !  coefficient PS         ! 
    !=========================!
    !$OMP do
    do me = 1, N_inn
        V_tmp(me) = sqrt(SP_uvw(me,1)**2.0d0 + SP_uvw(me,2)**2.0d0)
    enddo
    !$OMP enddo
    !$OMP barrier

    !$OMP single
    V_max = maxval(V_tmp)
    coe_PS = PS_C * V_max * Delta_t_eff * h  ! coefficient for PS
    !$OMP end single
    !$OMP barrier

    !=========================!
    !  main                   ! 
    !=========================!
    !$OMP do
    do me = 1, N_inn

        !=========================!
        !  zero clear             ! 
        !=========================!
        N_link = 0
        x_link = 0.0d0
        y_link = 0.0d0
        r_link = 0.0d0
        W_link = 0.0d0
        dW_link  = 0.0d0
        rho_link = 0.0d0
        uvw_link = 0.0d0
        tem_link = 0.0d0
        d_x = 0.0d0
        d_r = 0.0d0
        d_u = 0.0d0
        d_v = 0.0d0
        d_T = 0.0d0

        my_cell = cell(me)

        do link_cell = 1, 9
            your_cell = cell_9(link_cell, my_cell)  ! neighboring 9 cells

            if ((your_cell < 1) .or. (cell_max < your_cell)) cycle ! out of range
            if (N_cell(your_cell) == 0) cycle                      ! no particle in your_cell

            do num = 1, N_cell(your_cell)
                you = cell_info(num, your_cell) ! your name

                if (me == you) cycle            ! if me=you, skip

                N_link = N_link + 1             ! number of linked neighbors

                if (N_link > N_you) stop '[error] N_link > N_you @ cal_PS.f90'

                x_link(N_link) = SP_xyz(you,1) - SP_xyz(me,1) ! x_ji = x_j - x_i
                y_link(N_link) = SP_xyz(you,2) - SP_xyz(me,2) ! y_ji = y_j - y_i
                r_link(N_link) = sqrt(x_link(N_link)**2.0d0 + y_link(N_link)**2.0d0)

                call cal_W (r_link(N_link), h,  W_link(N_link)) ! from lib_kernel_function
                call cal_dW(r_link(N_link), h, dW_link(N_link)) ! from lib_kernel_function
                rho_link(  N_link) = SP_rho(you)
                tem_link(  N_link) = SP_tem(you  ) - SP_tem(me  )
                uvw_link(:,N_link) = SP_uvw(you,:) - SP_uvw(me,:)
            enddo
        enddo

        !===================================!
        !  shift vector & Taylor expansion  ! 
        !===================================!
#if defined(PS_1ST)
        call cal_PS_vec(2, N_you, N_link, W_link(:), dW_link(:), SP_rho(me), &
                        x_link(:), y_link(:), r_link(:), rho_link(:), uvw_link(:,:), tem_link(:), &
                        d_x(:), d_r(1:2), d_u(1:2), d_v(1:2), d_T(1:2))

        d_x(:) = - coe_PS * d_x(:)

        if (sqrt(d_x(1)**2.0d0 + d_x(2)**2.0d0) > 0.2d0*h) then  ! for safety
            d_x(:) = 0.2d0*h * d_x(:) / sqrt(d_x(1)**2.0d0 + d_x(2)**2.0d0)
        endif

        SP_xyz_new(me,:) = SP_xyz(me,:) + d_x(:)

        SP_uvw_new(me,1) = SP_uvw(me,1) + (d_u(1)*d_x(1) + d_u(2)*d_x(2))
        SP_uvw_new(me,2) = SP_uvw(me,2) + (d_v(1)*d_x(1) + d_v(2)*d_x(2))
        SP_rho_new(me)   = SP_rho(me)   + (d_r(1)*d_x(1) + d_r(2)*d_x(2))
        SP_tem_new(me)   = SP_tem(me)   + (d_T(1)*d_x(1) + d_T(2)*d_x(2))

#elif defined(PS_2ND)
        call cal_PS_vec(5, N_you, N_link, W_link(:), dW_link(:), SP_rho(me), &
                        x_link(:), y_link(:), r_link(:), rho_link(:), uvw_link(:,:), tem_link(:), &
                        d_x(:), d_r(1:5), d_u(1:5), d_v(1:5), d_T(1:5))

        d_x(:) = - coe_PS * d_x(:)

        if (sqrt(d_x(1)**2.0d0 + d_x(2)**2.0d0) > 0.2d0*h) then  ! for safety
            d_x(:) = 0.2d0*h * d_x(:) / sqrt(d_x(1)**2.0d0 + d_x(2)**2.0d0)
        endif
        
        SP_xyz_new(me,:) = SP_xyz(me,:) + d_x(:)
        
        SP_uvw_new(me,1) = SP_uvw(me,1) &
                         + (d_u(1)*d_x(1) + d_u(2)*d_x(2)) &
                         + (d_u(3)*d_x(1)**2.0d0 + d_u(4)*d_x(1)*d_x(2)*2.0d0 + d_u(5)*d_x(2)**2.0d0) / 2.0d0

        SP_uvw_new(me,2) = SP_uvw(me,2) &
                         + (d_v(1)*d_x(1) + d_v(2)*d_x(2)) &
                         + (d_v(3)*d_x(1)**2.0d0 + d_v(4)*d_x(1)*d_x(2)*2.0d0 + d_v(5)*d_x(2)**2.0d0) / 2.0d0

        SP_rho_new(me)   = SP_rho(me) &
                         + (d_r(1)*d_x(1) + d_r(2)*d_x(2)) &
                         + (d_r(3)*d_x(1)**2.0d0 + d_r(4)*d_x(1)*d_x(2)*2.0d0 + d_r(5)*d_x(2)**2.0d0) / 2.0d0
        
        SP_tem_new(me)   = SP_tem(me)  &
                         + (d_T(1)*d_x(1) + d_T(2)*d_x(2)) &
                         + (d_T(3)*d_x(1)**2.0d0 + d_T(4)*d_x(1)*d_x(2)*2.0d0 + d_T(5)*d_x(2)**2.0d0) / 2.0d0
#endif

    enddo
    !$OMP enddo
    !$OMP barrier
    !$OMP end parallel

    SP_xyz(1:N_inn,:) = SP_xyz_new(:,:)
    SP_uvw(1:N_inn,:) = SP_uvw_new(:,:)
    SP_rho(1:N_inn)   = SP_rho_new(:)
    SP_tem(1:N_inn)   = SP_tem_new(:)

    call cal_outside                                      ! if out of range
    call cal_EOS(SP_rho(1:N_inn), SP_pre(1:N_inn), N_inn) ! cal pressure from EOS
    call cal_cell_info                                    ! cal cell_info

end subroutine cal_PS

! END !