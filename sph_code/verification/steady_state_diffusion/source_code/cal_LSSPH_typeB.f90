subroutine cal_LSSPH_typeB(q, N_you, N_link, m, h, x, W, r, f, ans)
    !=========================!
    !  local variables        !
    !=========================!
    implicit none
    integer, intent(in) :: q, N_you, N_link
    real(8), intent(in) :: m, h, x(N_you,3), W(N_you), r(N_you), f(N_you)
    real(8), intent(out) :: ans

    integer :: you, IPIV(q), INFO, i, j
    real(8) :: x_ji, y_ji, V, a_vec(q), M_mat(q,q), M_mat_cp(q,q), b_vec(q)
    
    a_vec = 0.0d0  ! zero clear
    M_mat = 0.0d0
    b_vec = 0.0d0

    do you = 1, N_link
        x_ji = - x(you,1)  ! x = x_i - x_j
        y_ji = - x(you,2)
        V    = m / r(you)  ! volume

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
                M_mat(i,j) = M_mat(i,j) + V * W(you) * a_vec(i) * a_vec(j)
            enddo
            b_vec(i) = b_vec(i) + V * W(you) * a_vec(i) * f(you)
        enddo
    enddo

    !=========================!
    !  solve Md=b             !
    !=========================!
    M_mat_cp(:,:) = M_mat(:,:)
    call DGESV(q, 1, M_mat_cp, q, IPIV, b_vec, q, INFO)
    ans = b_vec(1)

end subroutine cal_LSSPH_typeB

! END !