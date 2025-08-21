subroutine cal_VM_to_WL
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
    integer :: me, you, my_cell, your_cell, link_cell, num, N_link
    real(8) :: x_link(N_you), y_link(N_you), r_link(N_you), W_link(N_you)
    real(8) :: rho_link(N_you), pre_link(N_you), uvw_link(2,N_you), tem_link(N_you)

    !=========================!
    !  zero clear             ! 
    !=========================!
    VM_rho = 0.0d0
    VM_pre = 0.0d0
    VM_uvw = 0.0d0
    VM_tem = 0.0d0

    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me, you, my_cell, your_cell, link_cell, num, N_link) &
    !$OMP private(x_link, y_link, r_link, W_link, rho_link, pre_link, uvw_link, tem_link) &
    !$OMP shared(N_VM, cell_VM, cell_max, cell_9, N_cell, cell_info, h, N_you, temp_ave) &
    !$OMP shared(SP_kind, SP_xyz, SP_rho, SP_pre, SP_uvw, SP_tem) &
    !$OMP shared(VM_kind, VM_xyz, VM_rho, VM_pre, VM_uvw, VM_tem, VM_WL_name)

    !=========================!
    !  cal VM                 !
    !=========================!
    !$OMP do
    do me = 1, N_VM ! loop for Virtual Markers
        
        !=========================!
        !  zero clear             ! 
        !=========================!
        N_link = 0
        x_link = 0.0d0
        y_link = 0.0d0
        r_link = 0.0d0
        W_link = 0.0d0
        rho_link = 0.0d0
        pre_link = 0.0d0
        uvw_link = 0.0d0
        tem_link = 0.0d0

        my_cell = cell_VM(me) ! my cell number

        !=========================!
        !  cal link list          ! 
        !=========================!
        do link_cell = 1, 9
            your_cell = cell_9(link_cell, my_cell) ! neighboring 9 cells

            if ((your_cell < 1) .or. (cell_max < your_cell)) cycle ! out of range
            if (N_cell(your_cell) == 0) cycle                      ! no particles in your_cell

            do num = 1, N_cell(your_cell)
                you = cell_info(num, your_cell) ! your name
                
                if (SP_kind(you) == 0) then     ! if you are an inner particle
                    N_link = N_link + 1         ! number of linked neighbors
                    
                    if (N_link > N_you) stop '[error] N_link > N_you @ cal_VM_to_WL'

                    x_link(N_link) = VM_xyz(me,1) - SP_xyz(you,1) ! x_ij 
                    y_link(N_link) = VM_xyz(me,2) - SP_xyz(you,2) ! y_ij 
                    r_link(N_link) = sqrt(x_link(N_link)**2.0d0 + y_link(N_link)**2.0d0)
                    call cal_W(r_link(N_link), h, W_link(N_link)) ! from lib_kernel_function
                    rho_link(  N_link) = SP_rho(you)
                    pre_link(  N_link) = SP_pre(you)
                    uvw_link(:,N_link) = SP_uvw(you,:)
                    tem_link(  N_link) = SP_tem(you)
                endif
            enddo
        enddo

        !=========================!
        !  SPH approx for VM      ! 
        !=========================!
        if (N_link == 0) cycle ! if there is no neighbors, skip

#if defined(WL_0TH)
        call cal_WL_0th(N_you, N_link, W_link(:), &
                        rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), &
                        VM_rho(me), VM_pre(me), VM_uvw(me,1), VM_uvw(me,2), VM_tem(me))
#elif defined(WL_1ST)
        call cal_WL_1st(N_you, N_link, W_link(:), &
                        rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), &
                        VM_rho(me), VM_pre(me), VM_uvw(me,1), VM_uvw(me,2), VM_tem(me))
#elif defined(WL_2ND)
        call cal_WL_LSSPH(3, N_you, N_link, W_link(:), x_link(:), y_link(:), &
                          rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), &
                          VM_rho(me), VM_pre(me), VM_uvw(me,1), VM_uvw(me,2), VM_tem(me))
#elif defined(WL_3RD)
        call cal_WL_LSSPH(6, N_you, N_link, W_link(:), x_link(:), y_link(:), &
                          rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), &
                          VM_rho(me), VM_pre(me), VM_uvw(me,1), VM_uvw(me,2), VM_tem(me))
#elif defined(WL_4TH)
        call cal_WL_LSSPH(10, N_you, N_link, W_link(:), x_link(:), y_link(:), &
                          rho_link(:), pre_link(:), uvw_link(:,:), tem_link(:), &
                          VM_rho(me), VM_pre(me), VM_uvw(me,1), VM_uvw(me,2), VM_tem(me))
#endif
        
        !=========================!
        !  VM to WL               ! 
        !=========================!
        you = VM_WL_name(me)       ! pair name
        ! SP_rho(you) = VM_rho(me) ! mirroring density
        SP_rho(you) = rho_ref      ! constant density
        SP_pre(you) = VM_pre(me)   ! mirroring pressure


        ! bottom wall
        if (VM_kind(me)==1) then
#if defined(NO_SLIP_BOTTOM)
            SP_uvw(you,:) = - VM_uvw(me,:)
#elif defined(FREE_SLIP_BOTTOM)
            SP_uvw(you,1) = + VM_uvw(me,1)
            SP_uvw(you,2) = - VM_uvw(me,2)
#endif

#if defined(BOUSSINESQ_CONV)
            SP_tem(you) = 2.0d0 * T_bot - VM_tem(me)
            SP_pre(you) = SP_pre(you) + rho_ref* alpha* (T_bot-temp_ave)* gravity* &
                                        (SP_xyz(you,2) - VM_xyz(me,2))
#endif
        

        ! bottom corners wall
        elseif ((VM_kind(me)==5).or.(VM_kind(me)==6)) then
            SP_uvw(you,:) = - VM_uvw(me,:)
#if defined(BOUSSINESQ_CONV)
            SP_tem(you) = 2.0d0 * T_bot - VM_tem(me)
            SP_pre(you) = SP_pre(you) + rho_ref* alpha* (T_bot-temp_ave)* gravity* &
                                        (SP_xyz(you,2) - VM_xyz(me,2))
#endif


        ! top wall
        elseif (VM_kind(me)==2) then
#if defined(NO_SLIP_TOP)
            SP_uvw(you,:) = - VM_uvw(me,:)
#elif defined(FREE_SLIP_TOP)
            SP_uvw(you,1) = + VM_uvw(me,1)
            SP_uvw(you,2) = - VM_uvw(me,2)
#elif defined(CONSTANT_FLOW_TOP)
            SP_uvw(you,1) = 2.0d0*U_top - VM_uvw(me,1)
            SP_uvw(you,2) = - VM_uvw(me,2)
#endif

#if defined(BOUSSINESQ_CONV)
            SP_tem(you) = 2.0d0 * T_top - VM_tem(me)
            SP_pre(you) = SP_pre(you) + rho_ref* alpha* (T_top-temp_ave)* gravity* &
                                        (SP_xyz(you,2) - VM_xyz(me,2))
#endif


        ! top corners wall
        elseif ((VM_kind(me)==7).or.(VM_kind(me)==8)) then
#if defined(CONSTANT_FLOW_TOP)
            SP_uvw(you,1) = 2.0d0*U_top - VM_uvw(me,1)
            SP_uvw(you,2) = - VM_uvw(me,2)
#else
            SP_uvw(you,:) = - VM_uvw(me,:)
#endif

#if defined(BOUSSINESQ_CONV)
            SP_tem(you) = 2.0d0 * T_top - VM_tem(me)
            SP_pre(you) = SP_pre(you) + rho_ref* alpha* (T_top-temp_ave)* gravity* &
                                        (SP_xyz(you,2) - VM_xyz(me,2))
#endif

        
        ! left wall
        elseif (VM_kind(me)==3) then
#if defined(NO_SLIP_LEFT)
            SP_uvw(you,:) = - VM_uvw(me,:)
#elif defined(FREE_SLIP_LEFT)
            SP_uvw(you,1) = - VM_uvw(me,1)
            SP_uvw(you,2) = + VM_uvw(me,2)
#endif
#if defined(BOUSSINESQ_CONV)
            SP_tem(you) = VM_tem(me)
#endif


        ! right wall
        elseif (VM_kind(me)==4) then
#if defined(NO_SLIP_RIGHT)
            SP_uvw(you,:) = - VM_uvw(me,:)
#elif defined(FREE_SLIP_RIGHT)
            SP_uvw(you,1) = - VM_uvw(me,1)
            SP_uvw(you,2) = + VM_uvw(me,2)
#endif
#if defined(BOUSSINESQ_CONV)
            SP_tem(you) = VM_tem(me)
#endif

        else
            stop ' [error] VM_kind is not correct @ cal_VM_to_WL'
        endif
    enddo
    !$OMP enddo
    !$OMP end parallel

end subroutine cal_VM_to_WL

! END !