module lib_solve_equation
    contains
    !=========================!
    !  solve Ax=b             ! 
    !=========================!
    subroutine gauss_jordan_pv(A, x, b, n)
        implicit none
        integer, intent(in)  :: n                   ! size of A0(n,n), b(n)
        real(8), intent(in)  :: A(n,n), b(n)        ! linear equations Ax=b
        real(8), intent(out) :: x(n)                ! solution
        integer i, k, m
        real(8) ar, am, t, A_tmp(n,n), w(n)
        A_tmp(:,:) = A(:,:)                         ! copy origiranl A
        x(:)       = b(:)                           ! copy original b
        
        do k = 1, n
            !>>>>>>>>>>>>>>>>>>>>>
            ! partial pivoting
            !>>>>>>>>>>>>>>>>>>>>>
            m = k                                   ! row m including max pivot
            am = abs(A_tmp(k,k))                    ! original pivot
            do i = (k+1), n                         ! column i > (k+1) loop 
                if (abs(A_tmp(i,k)) > am) then
                    am = abs(A_tmp(i,k))            ! tmp max pivot
                    m  = i                          ! tmp row m including tmp max pivot 'am'
                endif
            enddo

            if (am == 0.0d0) stop 'error: A is singluar in gauss_jordan'

            ! change row k and row m
            if (k /= m) then
                w    (   k:n) = A_tmp(k, k:n)       ! escape row k
                A_tmp(k, k:n) = A_tmp(m, k:n)       ! substitute row m with max pivot into row k
                A_tmp(m, k:n) = w    (   k:n)       ! row k into row m
                t    = x(k)                         ! escape right-hand in row k
                x(k) = x(m)                         ! substitute row m with max pivot into row k
                x(m) = t                            ! row k into row m
            endif

            !>>>>>>>>>>>>>>>>>>>>>
            ! Gauss-Jordan
            !>>>>>>>>>>>>>>>>>>>>>
            ar = 1.0d0 / A_tmp(k,k)                 ! inverse pivot
            A_tmp(k,k) = 1.0d0                      ! set 1 in diagonal component in row k
            A_tmp(k, k+1:n) = A_tmp(k, k+1:n) * ar  ! ar times other coefficients in row k
            x(k) = x(k) * ar                        ! ar times right-hand in row k
            
            ! calculation for each row i
            do i = 1, n
                if (i /= k) then
                    ! A_tmp(i,k)      = coefficient in the same column k including pivot
                    ! A_tmp(k, k+1:n) = each coefficient in row k including pivot
                    A_tmp(i, k+1:n) = A_tmp(i, k+1:n) - A_tmp(i,k) * A_tmp(k, k+1:n)
                    x(i)            = x(i)            - A_tmp(i,k) * x(k)
                    A_tmp(i, k)     = 0.0d0
                endif
            enddo 
        enddo
    end subroutine gauss_jordan_pv

    !=========================!
    !  solve A^-1             ! 
    !=========================!
    subroutine gauss_jordan_inv(A, A_inv, n)
        implicit none
        integer, intent(in)  :: n                  ! size of A0(n,n)
        real(8), intent(in)  :: A(n,n)             ! linear equations Ax=b
        real(8), intent(out) :: A_inv(n,n)         ! solution A^(-1)

        integer i, k, m
        real(8) ar, am, A_tmp(n, 2*n), w(2*n)   ! A_tmp: augmented matrix

        A_tmp(:,:) = 0.0d0
        A_tmp(:,1:n) = A(:,:)     ! copy origiranl A

        do i = 1, n
            A_tmp(i,n+i) = 1.0d0  ! identity matrix
        enddo

        do k = 1, n
            !>>>>>>>>>>>>>>>>>>>>>
            ! partial pivoting
            !>>>>>>>>>>>>>>>>>>>>>
            m = k                        ! row m including max pivot
            am = abs(A_tmp(k,k))         ! original pivot
            do i = (k+1), n              ! column i > (k+1) loop 
                if (abs(A_tmp(i,k)) > am) then
                    am = abs(A_tmp(i,k)) ! tmp max pivot
                    m  = i               ! tmp row m including tmp max pivot 'am'
                endif
            enddo

            if (am == 0.0d0) stop '[STOP] A is singluar in gauss_jordan_inv'

            ! change row k and row m
            if (k /= m) then
                w    (   k:2*n) = A_tmp(k, k:2*n) ! escape row k
                A_tmp(k, k:2*n) = A_tmp(m, k:2*n) ! substitute row m with max pivot into row k
                A_tmp(m, k:2*n) = w    (   k:2*n) ! row k into row m
            endif

            !>>>>>>>>>>>>>>>>>>>>>
            ! Gauss-Jordan
            !>>>>>>>>>>>>>>>>>>>>>
            ar = 1.0d0 / A_tmp(k,k)                     ! inverse pivot
            A_tmp(k,k) = 1.0d0                          ! set 1 in diagonal component in row k
            A_tmp(k, k+1:2*n) = A_tmp(k, k+1:2*n) * ar  ! ar times other coefficients in row k
            
            ! calculation for each row i
            do i = 1, n
                if (i /= k) then
                    ! A_tmp(i,k)        = coefficient in the same column k including pivot
                    ! A_tmp(k, k+1:2*n) = each coefficient in row k including pivot
                    A_tmp(i, k+1:2*n) = A_tmp(i, k+1:2*n) - A_tmp(i,k) * A_tmp(k, k+1:2*n)
                    A_tmp(i,k) = 0.0d0
                endif
            enddo 
        enddo

        A_inv(:,:) = A_tmp(:, n+1:2*n)

    end subroutine gauss_jordan_inv
end module lib_solve_equation

! END !