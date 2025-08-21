subroutine cal_classical(n, N_link, rho_me, pre_me, &
                         W_link, x_link, y_link, r_link, &
                         rho_link, pre_link, uvw_link, tem_link, &
                         d_r, d_p, d_u, d_v, d_T)

    !=========================!
    !  module                 ! 
    !=========================!
    use input
    use global_variables

    !=========================!
    !  local variables        ! 
    !=========================!
    implicit none
    integer, intent(in) :: n, N_link
    real(8), intent(in) :: rho_me, pre_me
    real(8), intent(in) :: W_link(n), x_link(n), y_link(n), r_link(n)
    real(8), intent(in) :: rho_link(n), pre_link(n), uvw_link(2,n), tem_link(n)
    ! 1=dx, 2=dy, 3=Laplacian (1=Laplacian for rho and tem)
    real(8), intent(out) :: d_r(1), d_p(2), d_u(3), d_v(3), d_T(1)

    integer :: you
    real(8) :: V_j, dW_dx(2)

    ! zero clear
    d_r = 0.0d0
    d_p = 0.0d0
    d_u = 0.0d0
    d_v = 0.0d0
    d_T = 0.0d0
    dW_dx = 0.0d0

    do you = 1, N_link
        V_j = SP_mass / rho_link(you)
        dW_dx(1) = x_link(you) / r_link(you) * W_link(you) ! dW/dx
        dW_dx(2) = y_link(you) / r_link(you) * W_link(you) ! dW/dy

        !=========================!
        !  d_r (Laplacian)        ! 
        !=========================!
        ! classical Laplacian model
        d_r(1) = d_r(1) + 2.0d0 * V_j * W_link(you)/r_link(you) * (rho_me-rho_link(you))

        !=========================!
        !  d_T (Laplacian)        ! 
        !=========================!
#if defined(BOUSSINESQ_CONV)
        ! classical Laplacian model
        d_T(1) = d_T(1) + 2.0d0 * V_j * W_link(you)/r_link(you) * (-tem_link(you))
#endif
        !=========================!
        !  du & dv                ! 
        !=========================!
        ! classical difference model
        d_u(1:2) = d_u(1:2) + V_j * dW_dx(:) * uvw_link(1,you)
        d_v(1:2) = d_v(1:2) + V_j * dW_dx(:) * uvw_link(2,you)

        ! classical Laplacian model
        d_u(3) = d_u(3) + 2.0d0 * V_j * W_link(you)/r_link(you) * (-uvw_link(1,you))
        d_v(3) = d_v(3) + 2.0d0 * V_j * W_link(you)/r_link(you) * (-uvw_link(2,you))

        !=========================!
        !  d_p                    ! 
        !=========================!
#if defined(CLASSICAL_SUM)
        ! classical sum model
        d_p(:) = d_p(:) + SP_mass * dW_dx(:) &
                    * (pre_link(you)/rho_link(you)**2.0d0 + pre_me/rho_me**2.0d0)

#elif defined(CLASSICAL_DIF)
        ! classical difference model
        d_p(:) = d_p(:) + V_j * dW_dx(:) * (pre_link(you)- pre_me)
#endif
    enddo

end subroutine cal_classical

! END !