subroutine cal_EOS(rho, pre, N_array)
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
    integer, intent(in) :: N_array
    real(8), intent(in) :: rho(N_array)
    real(8), intent(out) :: pre(N_array)

    integer :: me
    
    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me) &
    !$OMP shared(rho, pre, sound_speed, N_array)
    !$OMP do
    do me = 1, N_array
#if defined(CAVITY_FLOW) || defined(TAYLOR_GREEN)
       pre(me) = (sound_speed/zeta)**2.0d0 * (rho(me) - rho_ref) ! p=c^2*delta_rho
       
#elif defined(BOUSSINESQ_CONV)
       pre(me) = sound_speed**2.0d0 * (rho(me) - rho_ref)        ! p=c^2*delta_rho
#endif
    enddo
    !$OMP enddo
    !$OMP end parallel

end subroutine cal_EOS

! END !