subroutine cal_Runge_Kutta
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
    integer :: RK

    ! n step
    real(8), allocatable :: SP_xyz_old(:,:), SP_uvw_old(:,:), SP_rho_old(:), SP_tem_old(:)
    ! (n+1) step
    real(8), allocatable :: SP_xyz_new(:,:), SP_uvw_new(:,:), SP_rho_new(:), SP_tem_new(:)
    real(8), allocatable :: EOC(:), EOM(:,:), EOE(:)

    !=========================!
    !  allocate               ! 
    !=========================!
    allocate(SP_xyz_old(N_inn,2), SP_uvw_old(N_inn,2), SP_rho_old(N_inn), SP_tem_old(N_inn))
    allocate(SP_xyz_new(N_inn,2), SP_uvw_new(N_inn,2), SP_rho_new(N_inn), SP_tem_new(N_inn))
    allocate(EOC(N_inn), EOM(N_inn,2), EOE(N_inn))

    !=========================!
    !  copy original value    ! 
    !=========================!
    SP_xyz_old(:,:) = SP_xyz(1:N_inn,:) ! save n step values
    SP_uvw_old(:,:) = SP_uvw(1:N_inn,:)
    SP_rho_old(:)   = SP_rho(1:N_inn)
    SP_tem_old(:)   = SP_tem(1:N_inn)

    SP_xyz_new(:,:) = SP_xyz(1:N_inn,:) ! (n+1) step values
    SP_uvw_new(:,:) = SP_uvw(1:N_inn,:)
    SP_rho_new(:)   = SP_rho(1:N_inn)
    SP_tem_new(:)   = SP_tem(1:N_inn)

    !=========================!
    !  Runge-Kutta method     ! 
    !=========================!
    do RK = 1, 4

        !=========================!
        !  cal EOC & EOM          ! 
        !=========================!
        EOC = 0.0d0 ! zero clear
        EOM = 0.0d0
        EOE = 0.0d0

        call cal_EOC_EOM_EOE(EOC, EOM, EOE, N_inn)     ! see cal_EOC_EOM_EOE.f90

        !=========================!
        !  renew values           !
        !=========================!
        if (RK==1) then
            ! n + k1
            SP_xyz_new(:,:) = SP_xyz_new(:,:) + SP_uvw(1:N_inn,:)      *(Delta_t_eff/6.0d0)
            SP_uvw_new(:,:) = SP_uvw_new(:,:) + EOM(:,:)               *(Delta_t_eff/6.0d0)
            SP_rho_new(:)   = SP_rho_new(:)   + EOC(:)                 *(Delta_t_eff/6.0d0)
            SP_tem_new(:)   = SP_tem_new(:)   + EOE(:)                 *(Delta_t_eff/6.0d0)
            ! cal k2
            SP_xyz(1:N_inn,:) = SP_xyz_old(:,:) + SP_uvw(1:N_inn,:)    *(Delta_t_eff/2.0d0)
            SP_uvw(1:N_inn,:) = SP_uvw_old(:,:) + EOM(:,:)             *(Delta_t_eff/2.0d0)
            SP_rho(1:N_inn  ) = SP_rho_old(:)   + EOC(:)               *(Delta_t_eff/2.0d0)
            SP_tem(1:N_inn  ) = SP_tem_old(:)   + EOE(:)               *(Delta_t_eff/2.0d0)
        elseif (RK==2) then
            ! n + k1 + k2
            SP_xyz_new(:,:) = SP_xyz_new(:,:) + 2.0d0*SP_uvw(1:N_inn,:)*(Delta_t_eff/6.0d0)
            SP_uvw_new(:,:) = SP_uvw_new(:,:) + 2.0d0*EOM(:,:)         *(Delta_t_eff/6.0d0)
            SP_rho_new(:)   = SP_rho_new(:)   + 2.0d0*EOC(:)           *(Delta_t_eff/6.0d0)
            SP_tem_new(:)   = SP_tem_new(:)   + 2.0d0*EOE(:)           *(Delta_t_eff/6.0d0)
            ! cal k3
            SP_xyz(1:N_inn,:) = SP_xyz_old(:,:) + SP_uvw(1:N_inn,:)    *(Delta_t_eff/2.0d0)
            SP_uvw(1:N_inn,:) = SP_uvw_old(:,:) + EOM(:,:)             *(Delta_t_eff/2.0d0)
            SP_rho(1:N_inn  ) = SP_rho_old(:)   + EOC(:)               *(Delta_t_eff/2.0d0)
            SP_tem(1:N_inn  ) = SP_tem_old(:)   + EOE(:)               *(Delta_t_eff/2.0d0)
        elseif (RK==3) then
            ! n + k1 + k2 + k3
            SP_xyz_new(:,:) = SP_xyz_new(:,:) + 2.0d0*SP_uvw(1:N_inn,:)*(Delta_t_eff/6.0d0)
            SP_uvw_new(:,:) = SP_uvw_new(:,:) + 2.0d0*EOM(:,:)         *(Delta_t_eff/6.0d0)
            SP_rho_new(:)   = SP_rho_new(:)   + 2.0d0*EOC(:)           *(Delta_t_eff/6.0d0)
            SP_tem_new(:)   = SP_tem_new(:)   + 2.0d0*EOE(:)           *(Delta_t_eff/6.0d0)
            ! cal k4
            SP_xyz(1:N_inn,:) = SP_xyz_old(:,:) + SP_uvw(1:N_inn,:)    *(Delta_t_eff)
            SP_uvw(1:N_inn,:) = SP_uvw_old(:,:) + EOM(:,:)             *(Delta_t_eff)
            SP_rho(1:N_inn  ) = SP_rho_old(:)   + EOC(:)               *(Delta_t_eff)
            SP_tem(1:N_inn  ) = SP_tem_old(:)   + EOE(:)               *(Delta_t_eff)
        elseif (RK==4) then
            ! n + k1 + k2 + k3 + k4
            SP_xyz_new(:,:) = SP_xyz_new(:,:) + SP_uvw(1:N_inn,:)      *(Delta_t_eff/6.0d0)
            SP_uvw_new(:,:) = SP_uvw_new(:,:) + EOM(:,:)               *(Delta_t_eff/6.0d0)
            SP_rho_new(:)   = SP_rho_new(:)   + EOC(:)                 *(Delta_t_eff/6.0d0)
            SP_tem_new(:)   = SP_tem_new(:)   + EOE(:)                 *(Delta_t_eff/6.0d0)
            ! renewal
            SP_xyz(1:N_inn,:) = SP_xyz_new(:,:)
            SP_uvw(1:N_inn,:) = SP_uvw_new(:,:)
            SP_rho(1:N_inn)   = SP_rho_new(:)
            SP_tem(1:N_inn)   = SP_tem_new(:)
        endif


        !=========================!
        !  if out of range        !
        !=========================!
        call cal_outside

        !=========================!
        !  renew pressure by EOS  !
        !=========================!
        call cal_EOS(SP_rho(1:N_inn), SP_pre(1:N_inn), N_inn)

        !=========================!
        !  renew cell             !
        !=========================!
        call cal_cell_info

        !=========================!
        !  renew wall             !
        !=========================!
        call cal_VM_to_WL
    enddo

    !=========================!
    !  deallocate             ! 
    !=========================!
    deallocate(SP_xyz_old, SP_uvw_old, SP_rho_old, SP_tem_old)
    deallocate(SP_xyz_new, SP_uvw_new, SP_rho_new, SP_tem_new)
    deallocate(EOC, EOM, EOE)

end subroutine cal_Runge_Kutta

! END !