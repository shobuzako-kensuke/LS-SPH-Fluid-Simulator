subroutine check_options
!=========================!
!  SPH_model              !
!=========================!
#if defined(classical_SPH_Laplacian) || \
    defined(LS_SPH_2ND) || defined(LS_SPH_3RD)
#else
    stop ' [error] Compilation option is not correct. See SPH_model in Makefile.' 
#endif

!=========================!
!  kernel_function        !
!=========================!
#if defined(CUBIC_SPLINE) || defined(QUINTIC_SPLINE) || defined(WENDLAND_C2) || \
    defined(WENDLAND_C4)  || defined(WENDLAND_C6)
#else
    stop ' [error] Compilation option is not correct. See kernel_function in Makefile.' 
#endif

!=========================!
!  variable_density       !
!=========================!
#if defined(VARIABLE_DENSITY) || defined(CONSTANT_DENSITY)
#else
    stop ' [error] Compilation option is not correct. See variable_density in Makefile.' 
#endif

!=========================!
!  virtual_marker         !
!=========================!
#if defined(VM_ON) || defined(VM_OFF)
#else
    stop ' [error] Compilation option is not correct. See virtual_marker in Makefile.' 
#endif

!=========================!
!  wall_accuracy          !
!=========================!
#if defined(WL_1ST) || defined(WL_2ND) || defined(WL_3RD) || defined(WL_4TH)
#else
    stop ' [error] Compilation option is not correct. See wall_accuracy in Makefile.' 
#endif

!=========================!
!  wall_bottom            !
!=========================!
#if defined(DIRICHLET_BOTTOM) || defined(NEUMANN_BOTTOM)
#else
    stop ' [error] Compilation option is not correct. See wall_bottom in Makefile.' 
#endif

!=========================!
!  wall_top               !
!=========================!
#if defined(DIRICHLET_TOP) || defined(NEUMANN_TOP)
#else
    stop ' [error] Compilation option is not correct. See wall_top in Makefile.' 
#endif

!=========================!
!  wall_left              !
!=========================!
#if defined(DIRICHLET_LEFT) || defined(NEUMANN_LEFT)
#else
    stop ' [error] Compilation option is not correct. See wall_left in Makefile.' 
#endif

!=========================!
!  wall_right             !
!=========================!
#if defined(DIRICHLET_RIGHT) || defined(NEUMANN_RIGHT)
#else
    stop ' [error] Compilation option is not correct. See wall_right in Makefile.' 
#endif

end subroutine check_options

! END !