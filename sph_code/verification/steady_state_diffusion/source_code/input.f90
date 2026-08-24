module input
    implicit none
    !=========================!
    !  file name              !
    !=========================!
    character(len=999), parameter :: save_name = 'test'
    !=========================!
    !  General settings       !
    !=========================!
    integer, parameter :: N_threads  = 8         ! number of threads in OpenMP
    integer, parameter :: write_step = 100       ! writing interval
    !=========================!
    !  SPH settings           !
    !=========================!
    integer, parameter :: Nx = 10                ! number of particles along x axis
    real(8), parameter :: Nh = 1.2d0             ! h = Nh * Delta x
    !=========================!
    !  system parameters      !
    !=========================!
    real(8), parameter :: rho_ref   = 1.0d0      ! reference density [kg m-3]
    real(8), parameter :: x_rand    = 0.3d0      ! position perturbation
    real(8), parameter :: threshold = 1.0d-14    ! iteration threshold

end module input

! END !