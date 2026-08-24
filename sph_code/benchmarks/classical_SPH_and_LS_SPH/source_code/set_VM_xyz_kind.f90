subroutine set_VM_xyz_kind
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
    integer :: me, you
    
    !=========================!
    !  allocate               !
    !=========================!
    N_VM = N_out ! the number of Virtual Marker (VM) particles

    allocate(VM_kind(N_VM))    ! (name)
    allocate(VM_WL_name(N_VM)) ! VM_WL_name(me)=pair's name
    allocate(VM_xyz (N_VM, 2)) ! (name, (x,y))

    allocate(VM_rho(N_VM))     ! VM_rho(name)
    allocate(VM_pre(N_VM))     ! VM_pre(name)
    allocate(VM_uvw(N_VM, 2))  ! VM_uvw(name, (u,v))
    allocate(VM_tem(N_VM))     ! VM_tem(name)

    VM_kind = 0  ! zero clear
    VM_WL_name = 0
    VM_xyz = 0.0d0
    VM_rho = 0.0d0
    VM_pre = 0.0d0
    VM_uvw = 0.0d0
    VM_tem = 0.0d0

    !=========================!
    !  define position        !
    !=========================!
    me = 0 ! zero clear
    
    do you = N_inn+1, N_sys ! loop for wall particles

        me = me + 1          ! my name
        VM_WL_name(me) = you ! save your name

        if (SP_kind(you) == 1) then ! bottom wall
            VM_xyz(me,1) = SP_xyz(you,1)
            VM_xyz(me,2) = 2.0d0*(WL_thick + 0.0d0) - SP_xyz(you,2)
            VM_kind(me) = 1

        elseif (SP_kind(you) == 2) then  ! top wall
            VM_xyz(me,1) = SP_xyz(you,1)
            VM_xyz(me,2) = 2.0d0*(WL_thick + height) - SP_xyz(you,2)
            VM_kind(me) = 2

        elseif (SP_kind(you) == 3) then  ! left wall
            VM_xyz(me,1) = 2.0d0*(WL_thick + 0.0d0) - SP_xyz(you,1)
            VM_xyz(me,2) = SP_xyz(you,2)
            VM_kind(me) = 3

        elseif (SP_kind(you) == 4) then  ! right wall
            VM_xyz(me,1) = 2.0d0*(WL_thick + width) - SP_xyz(you,1)
            VM_xyz(me,2) = SP_xyz(you,2)
            VM_kind(me) = 4

        elseif (SP_kind(you) == 5) then  ! left bottom wall
            VM_xyz(me,1) = 2.0d0*(WL_thick + 0.0d0) - SP_xyz(you,1)
            VM_xyz(me,2) = 2.0d0*(WL_thick + 0.0d0) - SP_xyz(you,2)
            VM_kind(me) = 5

        elseif(SP_kind(you) == 6) then  ! right bottom wall
            VM_xyz(me,1) = 2.0d0*(WL_thick + width) - SP_xyz(you,1)
            VM_xyz(me,2) = 2.0d0*(WL_thick + 0.0d0) - SP_xyz(you,2)
            VM_kind(me) = 6

        elseif (SP_kind(you) == 7) then  ! left top wall
            VM_xyz(me,1) = 2.0d0*(WL_thick + 0.0d0)  - SP_xyz(you,1)
            VM_xyz(me,2) = 2.0d0*(WL_thick + height) - SP_xyz(you,2)
            VM_kind(me) = 7

        elseif (SP_kind(you) == 8) then  ! right top wall
            VM_xyz(me,1) = 2.0d0*(WL_thick + width)  - SP_xyz(you,1)
            VM_xyz(me,2) = 2.0d0*(WL_thick + height) - SP_xyz(you,2)
            VM_kind(me) = 8

        else
            stop '[error] inconsistent number of particles @ set_VN_xyz_kind'
        endif

    enddo

    ! for debug
    if (me /= N_out) stop '[error] N_VM is inconsistent @ set_VN_xyz_kind'
    
end subroutine set_VM_xyz_kind


! END !