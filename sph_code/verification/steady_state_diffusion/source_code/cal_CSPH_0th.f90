subroutine cal_CSPH_0th(N_you, N_link, m, W, r, f, ans)

    !=========================!
    !  local variables        !
    !=========================!
    implicit none
    integer, intent(in) :: N_you, N_link
    real(8), intent(in) :: m, W(N_you), r(N_you), f(N_you)
    real(8), intent(out) :: ans

    integer :: you
    real(8) :: V, sum, tmp
    sum = 0.0d0  ! zero clear
    tmp = 0.0d0

    !=========================!
    !  cal                    !
    !=========================!
    do you = 1, N_link
        V = m / r(you)                 ! volume
        tmp = tmp + W(you)* V* f(you)  ! function
        sum = sum + W(you)* V          ! renormalized term 
    enddo

    if (sum == 0.0d0) then
        sum = 1.0d0
    endif

    tmp = tmp / sum  ! renormalized by sum

    !=========================!
    !  substitution           !
    !=========================!
    ans = tmp

end subroutine cal_CSPH_0th

! END !