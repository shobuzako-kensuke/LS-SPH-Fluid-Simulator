subroutine cal_VM_to_WL
    !=========================!
    !  module                 !
    !=========================!
    !$use omp_lib
    use input
    use global_variables
    use lib_kernel_functions

    !=========================!
    !  local variables        !
    !=========================!
    implicit none
    integer :: me, you, my_cell, your_cell, link_cell, num, N_link
    real(8) :: pi=acos(-1.0d0), x(N_you,3), r(N_you), W(N_you), f(N_you)

    !=========================!
    !  zero clear             !
    !=========================!
    VM_f = 0.0d0

    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me, you, my_cell, your_cell, link_cell, num, N_link) &
    !$OMP private(x, r, W, f) &
    !$OMP shared(N_out, cell_VM, cell_max, cell_9, N_cell, cell_info, h, m, N_you) &
    !$OMP shared(SP_kind, SP_x, SP_r, SP_f) &
    !$OMP shared(VM_kind, VM_WL_name, VM_x, VM_f, pi, WL_thick)

    !=========================!
    !  cal VM                 !
    !=========================!
    !$OMP do
    do me = 1, N_out  ! loop for Virtual Markers
        
        !=========================!
        !  zero clear             !
        !=========================!
        N_link = 0
        x = 0.0d0
        r = 0.0d0
        W = 0.0d0
        f = 0.0d0

        my_cell = cell_VM(me)  ! my cell number

        !=========================!
        !  cal link list          !
        !=========================!
        do link_cell = 1, 9
            your_cell = cell_9(link_cell, my_cell)  ! neighboring 9 cells

            if ((your_cell < 1) .or. (cell_max < your_cell)) cycle  ! out of range
            if (N_cell(your_cell) == 0) cycle                       ! no particles in your_cell

            do num = 1, N_cell(your_cell)
                you = cell_info(num, your_cell)  ! your name
                
                if (SP_kind(you) == 0) then      ! if you are an inner particle
                    N_link = N_link + 1          ! number of linked neighbors
                    
                    if (N_link > N_you) stop ' [error] N_link > N_you @ cal_VM_to_WL.f90'

                    x(N_link,1:2) = VM_x(me,:) - SP_x(you,:)  ! x_ij = x_i - x_j
                    x(N_link,3  ) = sqrt(x(N_link,1)**2.0d0 + x(N_link,2)**2.0d0)
                    call cal_W(x(N_link,3), h, W(N_link))     ! from lib_kernel_functions
                    r(N_link) = SP_r(you)                     ! density (rho)
                    f(N_link) = SP_f(you)                     ! function value
                endif
            enddo
        enddo

        !=========================!
        !  SPH approx for VM      !
        !=========================!
        if (N_link == 0) cycle  ! if there is no neighbors, skip

#if defined(WL_1ST)
        call cal_CSPH_0th(N_you, N_link, m, W, r, f, VM_f(me))
#elif defined(WL_2ND)
        call cal_LSSPH_typeB(3, N_you, N_link, m, h, x, W, r, f, VM_f(me))
#elif defined(WL_3RD)
        call cal_LSSPH_typeB(6, N_you, N_link, m, h, x, W, r, f, VM_f(me))
#elif defined(WL_4TH)
        call cal_LSSPH_typeB(10, N_you, N_link, m, h, x, W, r, f, VM_f(me))
#endif
        
        !=========================!
        !  VM to WL               !
        !=========================!
        you = VM_WL_name(me)        ! pair name

        ! bottom wall
        if (VM_kind(me)==1) then

#if defined(DIRICHLET_BOTTOM)
            SP_f(you) = 2.0d0* (sin(pi* (VM_x(me,1)-WL_thick))) - VM_f(me)
#elif defined(NEUMANN_BOTTOM)
            SP_f(you) = VM_f(me) + 0.0d0* (SP_x(you,2) - VM_x(me,2))
#endif

        ! top wall
        elseif (VM_kind(me)==2) then
#if defined(DIRICHLET_TOP)
            SP_f(you) = 2.0d0* (-sin(pi* (VM_x(me,1)-WL_thick))) - VM_f(me)
#elif defined(NEUMANN_TOP)
            SP_f(you) = VM_f(me) + 0.0d0* (SP_x(you,2) - VM_x(me,2))
#endif

        ! left wall
        elseif (VM_kind(me)==3) then
#if defined(DIRICHLET_LEFT)
            SP_f(you) = 0.0d0 - VM_f(me)
#elif defined(NEUMANN_LEFT)
            SP_f(you) = VM_f(me) + pi*cos(pi* (VM_x(me,2)-WL_thick)) * (SP_x(you,1) - VM_x(me,1))
#endif

        ! right wall
        elseif (VM_kind(me)==4) then
#if defined(DIRICHLET_RIGHT)
            SP_f(you) = 0.0d0 - VM_f(me)
#elif defined(NEUMANN_RIGHT)
            SP_f(you) = VM_f(me) + (-pi)*cos(pi* (VM_x(me,2)-WL_thick)) * (SP_x(you,1) - VM_x(me,1))
#endif

        ! corners wall >> always Dirichlet conditions (f=0)
        elseif ((VM_kind(me)==5).or.(VM_kind(me)==6).or. &
                (VM_kind(me)==7).or.(VM_kind(me)==8)) then
                    SP_f(you) = 0.0d0 - VM_f(me)

        else
            stop ' [error] VM_kind is not correct @ cal_VM_to_WL.f90'
        endif
    enddo
    !$OMP enddo
    !$OMP end parallel

end subroutine cal_VM_to_WL

! END !