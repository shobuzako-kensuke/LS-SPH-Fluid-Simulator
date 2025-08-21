subroutine cal_classical_Laplacian(N_you, N_link, m, x, W, r, f, ans)
    !=========================!
    !  local variables        !
    !=========================!
    implicit none
    integer, intent(in) :: N_you, N_link
    real(8), intent(in) :: m, x(N_you,3), W(N_you), r(N_you), f(N_you)
    real(8), intent(out) :: ans

    integer :: you
    real(8) :: V, tmp
     
    tmp = 0.0d0  ! zero clear

    !=========================!
    !  cal Laplacian          !
    !=========================!
    do you = 1, N_link
        V = m / r(you)  ! volume
        tmp = tmp + 2.0d0 * V * W(you)/x(you,3) * (-f(you))
    enddo

    ans = tmp  ! substitution

end subroutine cal_classical_Laplacian

! END !