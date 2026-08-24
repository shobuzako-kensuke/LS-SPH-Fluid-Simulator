subroutine set_system_info
    !=========================!
    !  module                 ! 
    !=========================!
    use input
    use global_variables

    !=========================!
    !  SPH info               ! 
    !=========================!
    Delta_x = width / dble(N_x)             ! uniform particle spacing
    SP_mass = rho_ref * (Delta_x)**2.0d0    ! mass for 2-dimensions
    h = N_h * Delta_x                       ! smoothing length
#if defined(QUINTIC_SPLINE) 
    h_eff = 3.0d0 * h                       ! h_eff for quintic spline
#else
    h_eff = 2.0d0 * h                       ! h_eff for the others
#endif
    N_WL = ceiling(h_eff / Delta_x)         ! number of wall particles
    WL_thick = N_WL * Delta_x               ! length of wall

    total_step = end_step - start_step + 1  ! including initial condition

    !=========================!
    !  physical info          ! 
    !=========================!
    Re = U_top * width / (vis_ref / rho_ref)      ! Reynolds number

    sound_speed   = sqrt(K_ref / rho_ref)         ! sound of speed [m s-1]
    kinematic_vis = vis_ref / rho_ref             ! kinematic viscosity [m2 s-1]
    thermal_dif   = k_th / (rho_ref * c_p)        ! thermal diffusivity [m2 s-1]
    
    Delta_temp     = T_bot - T_Top                 ! temperature difference
    temp_ave       = 0.5d0* (T_bot + T_top)        ! average temperature
    Ra = alpha* Delta_temp* gravity* height**3.0d0 / &
         (kinematic_vis * thermal_dif)            ! Rayleigh number
    Pr = kinematic_vis / thermal_dif              ! Prandtl  number

    !=========================!
    !  time step              ! 
    !=========================!
    Delta_t_CFL = coe_CFL * (h / sound_speed)           ! time step for CFL
    Delta_t_vis = coe_vis * (h**2.0d0 / kinematic_vis)  ! time step for momentum
    Delta_t_CFL_relax = Delta_t_CFL * zeta * xi         ! weakly compressible approx
    Delta_t_vis_relax = Delta_t_vis * xi**2.0d0         ! weakly compressible approx

#if defined(BOUSSINESQ_CONV)
    Delta_t_th = coe_vis * (h**2.0d0 / thermal_dif)     ! time step for energy
#else
    Delta_t_th = Delta_t_CFL_relax* 10.0d0              ! pseudo time step for energy
#endif

    ! effective time step
    Delta_t_eff = min(Delta_t_CFL_relax, Delta_t_vis_relax, Delta_t_th)

    write(*,*) '+ ------------------------------------------------------------------------ +'
    write(*,*) '[message] time step condition below'
    write(*,*) '          1. Courant-Friedrichs-Lewy condition :', Delta_t_CFL
    write(*,*) '             >> relax                          :', Delta_t_CFL_relax
    write(*,*) '          2. momentum diffusion condition      :', Delta_t_vis
    write(*,*) '             >> relax                          :', Delta_t_vis_relax
    write(*,*) '          3. thermal  diffusion condition      :', Delta_t_th
    write(*,*) ''
    write(*,*) '          ==> effective time step              :', Delta_t_eff
    write(*,*) '+ ------------------------------------------------------------------------ +'
    write(*,*) ''
    write(*,*) ''

end subroutine set_system_info

! END !