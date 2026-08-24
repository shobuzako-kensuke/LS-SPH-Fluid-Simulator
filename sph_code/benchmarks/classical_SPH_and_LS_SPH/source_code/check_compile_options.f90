subroutine check_compile_options
!=========================!
!  problem_name           !
!=========================!
#if   defined(CAVITY_FLOW) || defined(TAYLOR_GREEN) || defined(BOUSSINESQ_CONV)
#else
    stop ' [error] compile option is not correct. see problem_name in Makefile.' 
#endif
    
!=========================!
!  kernel_function        !
!=========================!
#if   defined(CUBIC_SPLINE) || defined(QUINTIC_SPLINE) || defined(WENDLAND_C2)
#elif defined(WENDLAND_C4)  || defined(WENDLAND_C6)
#else
    stop ' [error] compile option is not correct. see kernel_function in Makefile.' 
#endif

!=========================!
!  SPH_model              !
!=========================!
#if   defined(CLASSICAL_SUM) || defined(CLASSICAL_DIF)  || defined(CSPH)
#elif defined(LSSPH_2ND)|| defined(LSSPH_3RD)
#else
    stop ' [error] compile option is not correct. see SPH_model in Makefile.'
#endif

!=========================!
!  particle_shifting      !
!=========================!
#if   defined(PS_OFF) || defined(PS_1ST) || defined(PS_2ND)
#else
    stop '[error] compile option is not correct. see particle_shifting in Makefile.'
#endif

!=========================!
!  wall_model             !
!=========================!
#if   defined(WL_0TH) || defined(WL_1ST) || defined(WL_2ND)
#elif defined(WL_3RD) || defined(WL_4TH)
#else
    stop ' [error] compile option is not correct. see wall_model in Makefile.'
#endif

!=========================!
!  wall_uvw_top           !
!=========================!
#if   defined(NO_SLIP_TOP) || defined(FREE_SLIP_TOP) || defined(CONSTANT_FLOW_TOP)
#else
    stop ' [error] compile option is not correct. see wall_uvw_top in Makefile.'
#endif

!=========================!
!  wall_uvw_bottom        !
!=========================!
#if   defined(NO_SLIP_BOTTOM) || defined(FREE_SLIP_BOTTOM)
#else
    stop ' [error] compile option is not correct. see wall_uvw_bottom in Makefile.'
#endif

!=========================!
!  wall_uvw_left          !
!=========================!
#if   defined(NO_SLIP_LEFT) || defined(FREE_SLIP_LEFT)
#else
    stop ' [error] compile option is not correct. see wall_uvw_left in Makefile.'
#endif

!=========================!
!  wall_uvw_right         !
!=========================!
#if   defined(NO_SLIP_RIGHT) || defined(FREE_SLIP_RIGHT)
#else
    stop ' [error] compile option is not correct. see wall_uvw_right in Makefile.'
#endif

end subroutine check_compile_options

! END !