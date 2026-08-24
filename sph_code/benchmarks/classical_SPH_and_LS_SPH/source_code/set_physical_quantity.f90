subroutine set_physical_quantity
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
    integer :: me
    real(8) :: pi=acos(-1.0d0), x, y
    
    !=========================!
    !  allocate               ! 
    !=========================!
    allocate(SP_rho(N_sys))    ! density (rho)
    allocate(SP_pre(N_sys))    ! pressure
    allocate(SP_uvw(N_sys,2))  ! velocity(name, (u,v))
    allocate(SP_tem(N_sys))    ! temperature

    !=========================!
    !  cal W_ave for PS       !
    !=========================!
    W_ave = 0.0d0
    call cal_W(Delta_x, h, W_ave)
    
    !=========================!
    !  initial state          ! 
    !=========================!
    SP_rho = rho_ref  ! set reference density
    SP_pre = 0.0d0    ! set static pressure
    SP_uvw = 0.0d0    ! set static velocity
    SP_tem = 0.0d0    ! set zero temperature

#if defined(TAYLOR_GREEN)
    do me = 1, N_inn
        x = pi* (SP_xyz(me,1) - WL_thick) / width
        y = pi* (SP_xyz(me,2) - WL_thick) / width
    
        SP_uvw(me,1) =  sin(x)* cos(y)
        SP_uvw(me,2) = -cos(x)* sin(y)

        SP_pre(me) = 0.25d0* (cos(2.0d0*x) + cos(2.0d0*y))
        
        SP_rho(me) = rho_ref + SP_pre(me) / (sound_speed/zeta)**2.0d0
    enddo
#endif

#if defined(BOUSSINESQ_CONV)
    do me = 1, N_inn
        x = SP_xyz(me,1) - WL_thick  ! translation
        y = SP_xyz(me,2) - WL_thick  ! translation
        
        SP_tem(me) = - (Delta_temp / height)* y + T_bot  ! linear profile of temperature
        if (y <= 0.1d0* height) then
            SP_tem(me) = SP_tem(me) + 1.0d-2* Delta_temp* cos(pi* x/ width) ! perturbation
        endif
        ! SP_pre(me) = rho_ref* alpha* (SP_tem(me) - temp_ave)* gravity* y
        ! SP_rho(me) = rho_ref + SP_pre(me) / sound_speed**2.0d0
    enddo
#endif

end subroutine set_physical_quantity

! END !