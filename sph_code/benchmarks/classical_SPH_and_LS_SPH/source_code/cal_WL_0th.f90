subroutine cal_WL_0th(n, N_link, W_link, &
                      rho_link, pre_link, uvw_link, tem_link, &
                      ans_rho, ans_pre, ans_u, ans_v, ans_tem)

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
    real(8), intent(in) :: W_link(n), rho_link(n), pre_link(n), uvw_link(2,n), tem_link(n)
    real(8), intent(out) :: ans_rho, ans_pre, ans_u, ans_v, ans_tem

    integer :: you
    real(8) :: V_j, ans(5) ! rho, pre, u, v, tem
    ans = 0.0d0            ! zero clear

    !=========================!
    !  cal                    ! 
    !=========================!
    do you = 1, N_link
        V_j = SP_mass / rho_link(you)
        ans(1) = ans(1) + V_j * rho_link(  you) * W_link(you) ! rho
        ans(2) = ans(2) + V_j * pre_link(  you) * W_link(you) ! pre
        ans(3) = ans(3) + V_j * uvw_link(1,you) * W_link(you) ! u
        ans(4) = ans(4) + V_j * uvw_link(2,you) * W_link(you) ! v
#if defined(BOUSSINESQ_CONV)
        ans(5) = ans(5) + V_j * tem_link(  you) * W_link(you) ! tem
#endif
    enddo

    !=========================!
    !  substitution           ! 
    !=========================!
    ans_rho = ans(1)
    ans_pre = ans(2)
    ans_u   = ans(3)
    ans_v   = ans(4)
    ans_tem = ans(5)

end subroutine cal_WL_0th

! END !