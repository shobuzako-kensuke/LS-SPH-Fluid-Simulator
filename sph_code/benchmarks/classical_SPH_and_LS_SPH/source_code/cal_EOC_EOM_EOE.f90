subroutine cal_EOC_EOM_EOE(EOC, EOM, EOE, N_array)
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
    integer, intent(in) :: N_array
    real(8), intent(inout) :: EOC(N_array), EOM(N_array, 2), EOE(N_array)

    integer :: me, you, my_cell, your_cell, link_cell, num, N_link
    real(8) :: x_link(N_you), y_link(N_you), r_link(N_you), W_link(N_you)
    real(8) :: rho_link(N_you), pre_link(N_you), uvw_link(2,N_you), tem_link(N_you)
    real(8) :: d_r(5), d_p(5), d_u(5), d_v(5), d_T(5)

    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me, you, my_cell, your_cell, link_cell, num, N_link) &
    !$OMP private(x_link, y_link, r_link, W_link, rho_link, pre_link, uvw_link, tem_link) &
    !$OMP private(d_r, d_p, d_u, d_v, d_T) &
    !$OMP shared(N_inn, cell, cell_max, cell_9, N_cell, cell_info, h, N_you) &
    !$OMP shared(SP_kind, SP_xyz, SP_rho, SP_pre, SP_uvw, SP_tem) &
    !$OMP shared(sound_speed, temp_ave, kinematic_vis, thermal_dif, EOC, EOM, EOE)
    !$OMP do
    do me = 1, N_inn ! loop for inner particles

        !=========================!
        !  zero clear             ! 
        !=========================!
        N_link = 0
        x_link = 0.0d0
        y_link = 0.0d0
        r_link = 0.0d0
        W_link = 0.0d0
        rho_link = 0.0d0
        pre_link = 0.0d0
        uvw_link = 0.0d0
        tem_link = 0.0d0
        d_r = 0.0d0
        d_p = 0.0d0
        d_u = 0.0d0
        d_v = 0.0d0
        d_T = 0.0d0

        my_cell = cell(me)  ! my cell

        do link_cell = 1, 9
            your_cell = cell_9(link_cell, my_cell)  ! neighboring 9 cells

            if ((your_cell < 1) .or. (cell_max < your_cell)) cycle ! out of range
            if (N_cell(your_cell) == 0) cycle                      ! no particle in your_cell

            do num = 1, N_cell(your_cell)
                you = cell_info(num, your_cell) ! your name

                if (me == you) cycle ! if me=you, skip

                N_link = N_link + 1  ! number of linked neighbors

                if (N_link > N_you) stop '[error] N_link > N_you @ cal_EOC_EOM_EOE.f90'

                x_link(N_link) = SP_xyz(me,1) - SP_xyz(you,1) ! x_ij = x_i - x_j
                y_link(N_link) = SP_xyz(me,2) - SP_xyz(you,2) ! y_ij = y_i - y_j
                r_link(N_link) = sqrt(x_link(N_link)**2.0d0 + y_link(N_link)**2.0d0)

#if defined(LSSPH_2ND) || defined(LSSPH_3RD)
                call cal_W(r_link(N_link), h, W_link(N_link))  ! from lib_kernel_function
                rho_link(  N_link) = SP_rho(you)
                pre_link(  N_link) = SP_pre(you) - SP_pre(me)
                tem_link(  N_link) = SP_tem(you) - SP_tem(me)
                uvw_link(:,N_link) = SP_uvw(you,:) - SP_uvw(me,:)

#elif defined(CLASSICAL_SUM) || defined(CLASSICAL_DIF) || defined(CSPH)
                call cal_dW(r_link(N_link), h, W_link(N_link)) ! from lib_kernel_function
                rho_link(  N_link) = SP_rho(you)
                pre_link(  N_link) = SP_pre(you)
                tem_link(  N_link) = SP_tem(you  ) - SP_tem(me  )
                uvw_link(:,N_link) = SP_uvw(you,:) - SP_uvw(me,:)
#endif
            enddo
        enddo

        !=========================!
        !  SPH approximation      ! 
        !=========================!
        if (N_link == 0) cycle ! if there is no neighbors, skip

#if defined(CLASSICAL_SUM) || defined(CLASSICAL_DIF)
        call cal_classical(N_you, N_link, SP_rho(me), SP_pre(me), &
                           W_link(:), x_link(:), y_link(:), r_link(:), &
                           rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), &
                           d_r(1), d_p(1:2), d_u(1:3), d_v(1:3), d_T(1))
#elif defined(CSPH)
        call cal_CSPH(N_you, N_link, SP_rho(me), SP_pre(me), &
                      W_link(:), x_link(:), y_link(:), r_link(:), &
                      rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), &
                      d_r(1), d_p(1:2), d_u(1:3), d_v(1:3), d_T(1))
#elif defined(LSSPH_2ND)
        call cal_LSSPH(5, N_you, N_link, W_link(:), x_link(:), y_link(:), &
                       rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), SP_rho(me), &
                       d_r(:), d_p(:), d_u(:), d_v(:), d_T(:))
#elif defined(LSSPH_3RD)
        call cal_LSSPH(9, N_you, N_link, W_link(:), x_link(:), y_link(:), &
                       rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), SP_rho(me), &
                       d_r(:), d_p(:), d_u(:), d_v(:), d_T(:))
#endif

        !=========================!
        !  cal EOC, EOM, EOE      ! 
        !=========================!
#if defined(CAVITY_FLOW) || defined(TAYLOR_GREEN)
#if defined(CLASSICAL_SUM) || defined(CLASSICAL_DIF) || defined(CSPH)
        EOC(me  ) = -SP_rho(me)*(d_u(1)+d_v(2)) + coe_delta_SPH*(sound_speed/zeta)*h*d_r(1)
        EOM(me,1) = -d_p(1)/SP_rho(me) + (vis_ref/SP_rho(me))*d_u(3)
        EOM(me,2) = -d_p(2)/SP_rho(me) + (vis_ref/SP_rho(me))*d_v(3)
        
#elif defined(LSSPH_2ND) || defined(LSSPH_3RD)
        EOC(me  ) = -SP_rho(me)*(d_u(1)+d_v(2)) &
                  + coe_delta_SPH*(sound_speed/zeta)*h*(d_r(3)+d_r(5))
        EOM(me,1) = -d_p(1)/SP_rho(me) + (vis_ref/SP_rho(me))*(d_u(3)+d_u(5))
        EOM(me,2) = -d_p(2)/SP_rho(me) + (vis_ref/SP_rho(me))*(d_v(3)+d_v(5))
#endif


#elif defined(BOUSSINESQ_CONV)
#if defined(CLASSICAL_SUM) || defined(CLASSICAL_DIF) || defined(CSPH)
        EOC(me  ) = (- rho_ref * (d_u(1) + d_v(2))          ) / zeta**2.0d0
        EOM(me,1) = (-d_p(1)/rho_ref + kinematic_vis*d_u(3) ) / xi  **2.0d0  
        EOM(me,2) = (-d_p(2)/rho_ref + kinematic_vis*d_v(3) &
                     + alpha* (SP_tem(me)-temp_ave)* gravity) / xi  **2.0d0
        EOE(me  ) = thermal_dif * d_T(1)
        
#elif defined(LSSPH_2ND) || defined(LSSPH_3RD)
        EOC(me  ) = (- rho_ref * (d_u(1) + d_v(2))          ) / zeta**2.0d0
        EOM(me,1) = (- d_p(1)/rho_ref                       &
                     + kinematic_vis * (d_u(3) + d_u(5))    ) / xi  **2.0d0
        EOM(me,2) = (- d_p(2)/rho_ref                       &
                     + kinematic_vis * (d_v(3) + d_v(5))    &
                     + alpha* (SP_tem(me)-temp_ave)* gravity) / xi  **2.0d0
        EOE(me  ) = thermal_dif * (d_T(3) + d_T(5))
#endif
#endif

    enddo
    !$OMP enddo
    !$OMP end parallel

end subroutine cal_EOC_EOM_EOE

! END !