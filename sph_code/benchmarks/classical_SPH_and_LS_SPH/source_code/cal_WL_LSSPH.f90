subroutine cal_WL_LSSPH(q, n, N_link, W_link, x_link, y_link, &
                        rho_link, pre_link, uvw_link, tem_link, &
                        ans_rho, ans_pre, ans_u, ans_v, ans_tem)

    !=========================!
    !  module                 ! 
    !=========================!
    use input
    use global_variables
    use lib_solve_equation

    !=========================!
    !  local variables        ! 
    !=========================!
    implicit none
    integer, intent(in) :: q, n, N_link
    real(8), intent(in) :: x_link(n), y_link(n)
    real(8), intent(in) :: W_link(n), rho_link(n), pre_link(n), uvw_link(2,n), tem_link(n)
    real(8), intent(out) :: ans_rho, ans_pre, ans_u, ans_v, ans_tem

    integer :: you, IPIV(q), INFO, i, j
    real(8) :: x_ji, y_ji, V_j, a_vec(q), M_mat(q,q), M_mat_cp(q,q)
    real(8) :: b_vec(q,5), tmp_ans(q), ans(5) ! rho, pre, u, v, tem
    
    a_vec = 0.0d0 ! zero clear
    M_mat = 0.0d0
    b_vec = 0.0d0
    tmp_ans = 0.0d0
    ans     = 0.0d0

    do you = 1, N_link
        x_ji = - x_link(you)
        y_ji = - y_link(you)
        V_j = SP_mass / rho_link(you)

        !=========================!
        !  cal a_vec              ! 
        !=========================!
#if defined(WL_4TH)
        ! including 3rd-order derivatives
        a_vec(1) = 1.0d0
        a_vec(2) = x_ji                   / h
        a_vec(3) = y_ji                   / h
        a_vec(4) = x_ji**2.0d0     /2.0d0 / h**2.0d0
        a_vec(5) = x_ji*y_ji              / h**2.0d0
        a_vec(6) = y_ji**2.0d0     /2.0d0 / h**2.0d0
        a_vec(7) = x_ji**3.0d0     /6.0d0 / h**3.0d0
        a_vec(8) = x_ji**2.0d0*y_ji/2.0d0 / h**3.0d0
        a_vec(9) = x_ji*y_ji**2.0d0/2.0d0 / h**3.0d0
        a_vec(10)= y_ji**3.0d0     /6.0d0 / h**3.0d0

#elif defined(WL_3RD)
        ! including 2nd-order derivatives
        a_vec(1) = 1.0d0
        a_vec(2) = x_ji                   / h
        a_vec(3) = y_ji                   / h
        a_vec(4) = x_ji**2.0d0     /2.0d0 / h**2.0d0
        a_vec(5) = x_ji*y_ji              / h**2.0d0
        a_vec(6) = y_ji**2.0d0     /2.0d0 / h**2.0d0

#elif defined(WL_2ND)
        ! including 1st-order derivatives
        a_vec(1) = 1.0d0
        a_vec(2) = x_ji                   / h
        a_vec(3) = y_ji                   / h
#endif

        !=========================!
        !  cal M_mat & b_vec      ! 
        !=========================!
        do i = 1, q
            do j = 1, q
                M_mat(i,j) = M_mat(i,j) + V_j * W_link(you) * a_vec(i) * a_vec(j)
            enddo

            b_vec(i,1) = b_vec(i,1) + V_j * W_link(you) * a_vec(i) * rho_link(you)
            b_vec(i,2) = b_vec(i,2) + V_j * W_link(you) * a_vec(i) * pre_link(you)
            b_vec(i,3) = b_vec(i,3) + V_j * W_link(you) * a_vec(i) * uvw_link(1,you)
            b_vec(i,4) = b_vec(i,4) + V_j * W_link(you) * a_vec(i) * uvw_link(2,you)
#if defined(BOUSSINESQ_CONV)
            b_vec(i,5) = b_vec(i,5) + V_j * W_link(you) * a_vec(i) * tem_link(you)
#endif
        enddo
    enddo

    !=========================!
    !  solve Md=b             ! 
    !=========================!
#if defined(BOUSSINESQ_CONV)
    ! loop for (rho, pre, u, v, tem)
    do i = 1, 5
        tmp_ans (:)   = b_vec(:,i)
        M_mat_cp(:,:) = M_mat(:,:)
        call DGESV(q, 1, M_mat_cp, q, IPIV, tmp_ans, q, INFO)
        ans(i) = tmp_ans(1)
    enddo

#else
    ! loop for (rho, pre, u, v)
    do i = 1, 4
        tmp_ans (:)   = b_vec(:,i)
        M_mat_cp(:,:) = M_mat(:,:)
        call DGESV(q, 1, M_mat_cp, q, IPIV, tmp_ans, q, INFO)
        ans(i) = tmp_ans(1)
    enddo
#endif

    !=========================!
    !  substitution           ! 
    !=========================!
    ans_rho = ans(1)
    ans_pre = ans(2)
    ans_u   = ans(3)
    ans_v   = ans(4)
    ans_tem = ans(5)

end subroutine cal_WL_LSSPH

! END !