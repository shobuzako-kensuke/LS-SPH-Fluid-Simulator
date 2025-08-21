subroutine set_system
    !=========================!
    !  module                 !
    !=========================!
    !$use omp_lib
    use input
    use global_variables

    !=========================!
    !  local variables        !
    !=========================!
    implicit none
    real(8) :: rand_val

    integer :: me, Nx_fill, N_sys_fill
    real(8), allocatable :: x_sys(:,:), x_inn(:,:), x_out(:,:)

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                       fundamental information                            !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    D  = 1.0d0 / dble(Nx)                    ! Delta x [m]
    m  = 1.0d0* D**2.0d0                     ! mass = rho * (Delta x)**2 [kg]
    h  = Nh* D                               ! smoothing length [m]
#if defined(QUINTIC_SPLINE) 
    h_eff = 3.0d0 * h                        ! for quintic spline
#else
    h_eff = 2.0d0 * h                        ! otherwise
#endif

    N_WL  = ceiling(h_eff / D)               ! number along the wall thickness
    WL_thick = N_WL * D                      ! wall thickness [m]

    dt = 0.25d0* h**2.0d0 / 1.0d0             ! time step [s]    

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                        particle configuration                            !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !==============================!
    !  fill system with particles  !
    !==============================!
    Nx_fill    = Nx + 2*N_WL        ! inner + outer particles along x axis
    N_sys_fill = Nx_fill* Nx_fill   ! all particles if system packed

    allocate(x_sys(N_sys_fill, 2))  ! (name, (x,y))
    allocate(x_inn(N_sys_fill, 2))
    allocate(x_out(N_sys_fill, 2))

    x_sys = 0.0d0  ! zero clear
    x_inn = 0.0d0
    x_out = 0.0d0
    x_sys(1,:) = D / 2.0d0  ! first particle

    do me = 2, N_sys_fill
        if (me <= Nx_fill) then
            x_sys(me,1) = x_sys(me-1,1) + D
            x_sys(me,2) = x_sys(me-1,2)
        else
            x_sys(me,1) = x_sys(me-Nx_fill,1)
            x_sys(me,2) = x_sys(me-Nx_fill,2) + D
        endif
    enddo

    !=========================!
    !  categorize             !
    !=========================!
    N_inn = 0  ! zero clear
    N_out = 0

    do me = 1, N_sys_fill
        !=========================!
        !  find inner particles   !
        !=========================!
        if ((WL_thick + 0.0d0  < x_sys(me,1)) .and. &
            (WL_thick + 1.0d0  > x_sys(me,1)) .and. &
            (WL_thick + 0.0d0  < x_sys(me,2)) .and. &
            (WL_thick + 1.0d0  > x_sys(me,2))) then
            
            N_inn = N_inn + 1  ! count the number of inner particles
            x_inn(N_inn,:) = x_sys(me,:)

        !==============================!
        !  find wall particles         !
        !==============================!
        else
            N_out = N_out + 1  ! count the number of outer particles
            x_out(N_out,:) = x_sys(me,:)
        endif
    enddo

    !=========================!
    !  allocate & sort        !
    !=========================!
    N_sys = N_inn + N_out

    allocate(SP_kind(N_sys), SP_x(N_sys,2), SP_r(N_sys), SP_f(N_sys))
    allocate(SP_f_now(N_inn), SP_f_next(N_inn), RK(N_inn,2), dif_f(N_inn))
    SP_kind = 0  ! zero clear
    SP_x = 0.0d0
    SP_r = 0.0d0
    SP_f = 0.0d0
    SP_f_now  = 0.0d0
    SP_f_next = 0.0d0
    RK = 0.0d0
    dif_f = 0.0d0

    ! inner particles >> SP_kind = 0
    SP_x   (:N_inn,:) = x_inn(:N_inn,:)
    SP_kind(:N_inn  ) = 0

    ! wall particles >> SP_kind = 10 (tmp)
    SP_x   ((N_inn+1):,:) = x_out(:N_out,:)
    SP_kind((N_inn+1):  ) = 10
    
    !===========================!
    !  identify wall particles  !
    !===========================!
    do me = N_inn+1, N_sys

        ! bottom wall >> SP_kind = (1,5,6)
        if (SP_x(me,2) < WL_thick) then
            ! left bottom
            if (SP_x(me,1) < WL_thick) then
                SP_kind(me) = 5

            ! right bottom
            elseif (SP_x(me,1) > WL_thick + 1.0d0) then
                SP_kind(me) = 6
            
            ! bottom
            else
                SP_kind(me) = 1
            endif
        
        ! top wall >> SP_kind = (2,7,8)
        elseif (SP_x(me,2) > WL_thick + 1.0d0) then
            ! left top
            if (SP_x(me,1) < WL_thick) then
                SP_kind(me) = 7

            ! right top
            elseif (SP_x(me,1) > WL_thick + 1.0d0) then
                SP_kind(me) = 8
            
            ! top
            else
                SP_kind(me) = 2
            endif

        ! left wall >> SP_kind = 3
        elseif (SP_x(me,1) < WL_thick) then
            SP_kind(me) = 3

        ! right wall >> SP_kind = 4
        elseif (SP_x(me,1) > WL_thick + 1.0d0) then
            SP_kind(me) = 4

        else
            stop ' [error] SP_kind is not correct @ set_system.f90'
        endif
    enddo

    !=========================!
    !  deallocate             !
    !=========================!
    deallocate(x_sys, x_out, x_inn)

    !=========================!
    !  write                  !
    !=========================!
    write(*,*) '+ ------------------------------------------------------------------------ +'
    write(*,*) '[message] Nx, N_inn, N_out    :', Nx, N_inn, N_out
    write(*,*) '          Delta x         [m] :', D
    write(*,*) '          time step       [s] :', dt
    write(*,*) '          typical time step   :', 1.0d0 / dt
    write(*,*) '+ ------------------------------------------------------------------------ +'
    write(*,*) ''
    write(*,*) ''
    write(*,*) ''


    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                      virtual marker (VM) settings                        !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    call set_VM  ! set VM


    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                          position perturbation                           !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    call random_seed()
    do me = 1, N_inn
        call random_number(rand_val)
        SP_x(me,1) = SP_x(me,1)* (1.0d0 + x_rand* (rand_val - 0.5d0)* 2.0d0)
        call random_number(rand_val)
        SP_x(me,2) = SP_x(me,2)* (1.0d0 + x_rand* (rand_val - 0.5d0)* 2.0d0)
    enddo

end subroutine set_system

! END !