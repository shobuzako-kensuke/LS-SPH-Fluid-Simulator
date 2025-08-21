module input
    implicit none
    !=========================!
    !  file name              ! 
    !=========================!
    character(len=999), parameter :: save_name = 'TG_test'
    character(len=999), parameter :: read_name = 'new'
    !=========================!
    !  OpenMP                 ! 
    !=========================!
    integer, parameter :: N_threads = 8          ! number of threads in OpenMP
    !=========================!
    !  time integration       ! 
    !=========================!
    integer, parameter :: start_step = 1
    integer, parameter :: end_step   = 5000
    integer, parameter :: write_step = 100       ! writing interval
    real(8), parameter :: coe_CFL    = 0.5d0     ! Delta_t = coe_CFL * CFL
    real(8), parameter :: coe_vis    = 0.1d0     ! Delta_t = coe_vis * DIF
    real(8), parameter :: zeta       = 4.5d3     ! relaxation parameter for EOC
    real(8), parameter :: xi         = 1.0d0     ! relaxation parameter for EOM
    !=========================!
    !  SPH settings           ! 
    !=========================!
    integer, parameter :: N_x = 50               ! number of particles along x-axis
    integer, parameter :: N_y = 50               ! number of particles along y-axis
    real(8), parameter :: N_h = 1.2d0            ! h = N_h * Delta_x
    !=========================!
    !  physical parameter     ! 
    !=========================!
    real(8), parameter :: rho_ref = 1.0d0        ! reference density [kg m-3]
    real(8), parameter :: vis_ref = 1.0d-2       ! reference viscosity [Pa s]
    real(8), parameter :: K_ref   = 2.2d9        ! reference artificial bulk modulus [Pa]
    real(8), parameter :: k_th    = 5.0d0        ! thermal conductivity [W m-1 K-1]
    real(8), parameter :: c_p     = 1.25d3       ! specific heat [J kg-1 K-1]
    real(8), parameter :: alpha   = 2.5d-5       ! thermal expansion [K-1]
    !=========================!
    !  Particle Shifting      ! 
    !=========================!
    real(8), parameter :: PS_C = 1.5d0           ! dr=coe_PS*(1+PS_R*(W/W_ave)**PS_n)*(dW_ij)*V_j
    real(8), parameter :: PS_R = 0.2d0           ! coe_PS = PS_C * V_max * Delta t * h
    real(8), parameter :: PS_n = 4.0d0
    !=========================!
    !  system parameter       ! 
    !=========================!
    real(8), parameter :: width   = 1.0d0        ! system width  (along x-axis) [m]
    real(8), parameter :: height  = 1.0d0        ! system height (along y-axis) [m]
    real(8), parameter :: U_top   = 1.0d0        ! used for cavity flow
    real(8), parameter :: T_top   = 0.0d0        ! used for Boussinesq convection
    real(8), parameter :: T_bot   = 1.0d3
    real(8), parameter :: gravity = 1.0d1
    !=========================!
    !  delta-SPH              ! 
    !=========================!
    real(8), parameter :: coe_delta_SPH = 0.0d0  ! coefficient of density diffusion term

end module input

! END !