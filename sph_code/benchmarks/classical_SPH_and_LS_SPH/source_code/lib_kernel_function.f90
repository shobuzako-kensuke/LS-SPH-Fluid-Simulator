module lib_kernel_function
    implicit none
    contains
    !=========================!
    !  W                      ! 
    !=========================!
    subroutine cal_W(r, h, W)
        implicit none
        real(8), intent(in)  :: r, h
        real(8), intent(out) :: W
        real(8) :: q, pi = dacos(-1.0d0)
        q = r/h

#if defined (CUBIC_SPLINE)
        if ((0.0d0 <= q) .and. (q <= 1.0d0)) then
            W = (2.0d0-q)**3.0d0 - 4.0d0*(1.0d0-q)**3.0d0
        elseif ((1.0d0 < q) .and. (q <= 2.0d0)) then
            W = (2.0d0-q)**3.0d0
        else
            W = 0.0d0
        endif
        W = 5.0d0 / (14.0d0 * pi * h**2.0) * W
#endif

#if defined (QUINTIC_SPLINE)
        if ((0.0d0 <= q) .and. (q <= 1.0d0)) then
            W = (3.0d0-q)**5.0d0 - 6.0d0*(2.0d0-q)**5.0d0 + 15.0d0*(1.0d0-q)**5.0d0
        elseif ((1.0d0 < q) .and. (q <= 2.0d0)) then
            W = (3.0d0-q)**5.0d0 - 6.0d0*(2.0d0-q)**5.0d0
        elseif ((2.0d0 < q) .and. (q <= 3.0d0)) then
            W = (3.0d0-q)**5.0d0
        else
            W = 0.0d0
        endif
        W = 7.0d0 / (478.0d0 * pi * h**2.0d0) * W
#endif

#if defined (WENDLAND_C2)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            W = (1.0d0 - 0.5d0*q)**4.0d0 * (1.0d0 + 2.0d0*q)
        else
            W = 0.0d0
        endif
        W = 7.0d0 / (4.0d0 * pi * h**2.0d0) * W
#endif

#if defined (WENDLAND_C4)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            W = (1.0d0 - 0.5d0*q)**6.0d0 * (1.0d0 + 3.0d0*q + 35.0d0*q**2.0d0/12.0d0)
        else
            W = 0.0d0
        endif
        W = 9.0d0 / (4.0d0 * pi * h**2.0d0) * W
#endif

#if defined (WENDLAND_C6)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            W = (1.0d0 - 0.5d0*q)**8.0d0 * &
                (1.0d0 + 4.0d0*q + 25.0d0*q**2.0d0/4.0d0 + 4.0d0*q**3.0d0)
        else
            W = 0.0d0
        endif
        W = 39.0d0 / (14.0d0 * pi * h**2.0d0) * W
#endif
    end subroutine cal_W

    !=========================!
    !  dW                     ! 
    !=========================!
    subroutine cal_dW(r, h, dW)
        implicit none
        real(8), intent(in)  :: r, h
        real(8), intent(out) :: dW
        real(8) :: q, pi = dacos(-1.0d0)
        q = r/h

#if defined (CUBIC_SPLINE)
        if ((0.0d0 <= q) .and. (q <= 1.0d0)) then
            dW = (2.0d0-q)**2.0d0 - 4.0d0*(1.0d0-q)**2.0d0
        elseif ((1.0d0 < q) .and. (q <= 2.0d0)) then
            dW = (2.0d0-q)**2.0d0
        else
            dW = 0.0d0
        endif
        dW = 5.0d0 / (14.0d0 * pi * h**2.0d0) * dW * (-3.0d0/h)
#endif

#if defined (QUINTIC_SPLINE)
        if ((0.0d0 <= q) .and. (q <= 1.0d0)) then
            dW = (3.0d0-q)**4.0d0 - 6.0d0*(2.0d0-q)**4.0d0 + 15.0d0*(1.0d0-q)**4.0d0
        elseif ((1.0d0 < q) .and. (q <= 2.0d0)) then
            dW = (3.0d0-q)**4.0d0 - 6.0d0*(2.0d0-q)**4.0d0
        elseif ((2.0d0 < q) .and. (q <= 3.0d0)) then
            dW = (3.0d0-q)**4.0d0
        else
            dW = 0.0d0
        endif
        dW = 7.0d0 / (478.0d0 * pi * h**2.0d0) * dW *(-5.0d0/h)
#endif

#if defined (WENDLAND_C2)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            dW = -2.0d0*(1.0d0 - 0.5d0*q)**3.0d0 * (1.0d0 + 2.0d0*q) &
                 +2.0d0*(1.0d0 - 0.5d0*q)**4.0d0
        else
            dW = 0.0d0
        endif
        dW = 7.0d0 / (4.0d0 * pi * h**2.0d0) * dW / h
#endif

#if defined (WENDLAND_C4)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            dW = -3.0d0*(1.0d0 - 0.5d0*q)**5.0d0 * (1.0d0 + 3.0d0*q + 35.0d0*q**2.0d0/12.0d0) &
                 +1.0d0*(1.0d0 - 0.5d0*q)**6.0d0 * (3.0d0 + 35.0d0*q/6.0d0)
        else
            dW = 0.0d0
        endif
        dW = 9.0d0 / (4.0d0 * pi * h**2.0d0) * dW / h
#endif


#if defined (WENDLAND_C6)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            dW = -4.0d0*(1.0d0 - 0.5d0*q)**7.0d0 &
                * (1.0d0 + 4.0d0*q + 25.0d0*q**2.0d0/4.0d0 + 4.0d0*q**3.0d0) &
                +1.0d0*(1.0d0 - 0.5d0*q)**8.0d0 &
                * (4.0d0 + 25.0d0*q/2.0d0 + 12.0d0*q**2.0d0)
        else
            dW = 0.0d0
        endif
        dW = 39.0d0 / (14.0d0 * pi * h**2.0d0) * dW / h
#endif
    end subroutine cal_dW

    !=========================!
    !  ddW                    ! 
    !=========================!
    subroutine cal_ddW(r, h, ddW)
        implicit none
        real(8), intent(in)  :: r, h
        real(8), intent(out) :: ddW
        real(8) :: q, pi = dacos(-1.0d0)
        q = r/h

#if defined (CUBIC_SPLINE)
        if ((0.0d0 <= q) .and. (q <= 1.0d0)) then
            ddW = (2.0d0-q) - 4.0d0*(1.0d0-q)
        elseif ((1.0d0 < q) .and. (q <= 2.0d0)) then
            ddW = (2.0d0-q)
        else
            ddW = 0.0d0
        endif
        ddW = 5.0d0 / (14.0d0 * pi * h**2.0d0) * ddW * (6.0d0/h**2.0d0)
#endif

#if defined (QUINTIC_SPLINE)
        if ((0.0d0 <= q) .and. (q <= 1.0d0)) then
            ddW = (3.0d0-q)**3.0d0 - 6.0d0*(2.0d0-q)**3.0d0 + 15.0d0*(1.0d0-q)**3.0d0
        elseif ((1.0d0 < q) .and. (q <= 2.0d0)) then
            ddW = (3.0d0-q)**3.0d0 - 6.0d0*(2.0d0-q)**3.0d0
        elseif ((2.0d0 < q) .and. (q <= 3.0d0)) then
            ddW = (3.0d0-q)**3.0d0
        else
            ddW = 0.0d0
        endif
        ddW = 7.0d0 / (478.0d0 * pi * h**2.0d0) * ddW *(20.0d0/h**2.0d0)
#endif

#if defined (WENDLAND_C2)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            ddW = +3.0d0*(1.0d0 - 0.5d0*q)**2.0d0 * (1.0d0 + 2.0d0*q) &
                  -4.0d0*(1.0d0 - 0.5d0*q)**3.0d0 &
                  -4.0d0*(1.0d0 - 0.5d0*q)**3.0d0
        else
            ddW = 0.0d0
        endif
        ddW = 7.0d0 / (4.0d0 * pi * h**2.0d0) * ddW / h**2.0d0
#endif

#if defined (WENDLAND_C4)
if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
        ddW = +15.0d0/2.0d0*(1.0d0 - 0.5d0*q)**4.0d0 &
                           * (1.0d0 + 3.0d0*q + 35.0d0*q**2.0d0/12.0d0) &
              -3.0d0       *(1.0d0 - 0.5d0*q)**5.0d0 * (3.0d0 + 35.0d0*q/6.0d0) &
              -3.0d0       *(1.0d0 - 0.5d0*q)**5.0d0 * (3.0d0 + 35.0d0*q/6.0d0) &
              +1.0d0       *(1.0d0 - 0.5d0*q)**6.0d0 * (35.0d0/6.0d0)
            else
            ddW = 0.0d0
        endif
        ddW = 9.0d0 / (4.0d0 * pi * h**2.0d0) * ddW / h**2.0d0
#endif

#if defined (WENDLAND_C6)
        if ((0.0d0 <= q) .and. (q <= 2.0d0)) then
            ddW = +14.0d0*(1.0d0 - 0.5d0*q)**6.0d0 &
                         * (1.0d0 + 4.0d0*q + 25.0d0*q**2.0d0/4.0d0 + 4.0d0*q**3.0d0) &
                  -4.0d0 *(1.0d0 - 0.5d0*q)**7.0d0 * (4.0d0 + 25.0d0*q/2.0d0 + 12.0d0*q**2.0d0) &
                  -4.0d0 *(1.0d0 - 0.5d0*q)**7.0d0 * (4.0d0 + 25.0d0*q/2.0d0 + 12.0d0*q**2.0d0) &
                  +1.0d0 *(1.0d0 - 0.5d0*q)**8.0d0 * (25.0d0/2.0d0 + 24.0d0*q)
        else
            ddW = 0.0d0
        endif
        ddW = 39.0d0 / (14.0d0 * pi * h**2.0d0) * ddW / h**2.0d0
#endif
    end subroutine cal_ddW
end module lib_kernel_function

! END !