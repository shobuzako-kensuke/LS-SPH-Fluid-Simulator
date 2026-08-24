subroutine set_VM
    !=========================!
    !  module                 !
    !=========================!
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
    allocate(VM_kind   (N_out))  ! (name)
    allocate(VM_WL_name(N_out))  ! VM_WL_name(me) = pair's name
    allocate(VM_x(N_out, 2))     ! (name, (x,y))
    allocate(VM_f(N_out))        ! (name)

    VM_kind = 0                  ! zero clear
    VM_WL_name = 0
    VM_x = 0.0d0
    VM_f = 0.0d0

    !=========================!
    !  define position        !
    !=========================!
    me = 0  ! zero clear
    
    do you = N_inn+1, N_sys  ! loop for wall particles

        me = me + 1           ! my name
        VM_WL_name(me) = you  ! save your name

        if (SP_kind(you) == 1) then  ! bottom wall
            VM_x(me,1) = SP_x(you,1)
            VM_x(me,2) = 2.0d0*(WL_thick + 0.0d0) - SP_x(you,2)
            VM_kind(me) = 1

        elseif (SP_kind(you) == 2) then  ! top wall
            VM_x(me,1) = SP_x(you,1)
            VM_x(me,2) = 2.0d0*(WL_thick + 1.0d0) - SP_x(you,2)
            VM_kind(me) = 2

        elseif (SP_kind(you) == 3) then  ! left wall
            VM_x(me,1) = 2.0d0*(WL_thick + 0.0d0) - SP_x(you,1)
            VM_x(me,2) = SP_x(you,2)
            VM_kind(me) = 3

        elseif (SP_kind(you) == 4) then  ! right wall
            VM_x(me,1) = 2.0d0*(WL_thick + 1.0d0) - SP_x(you,1)
            VM_x(me,2) = SP_x(you,2)
            VM_kind(me) = 4

        elseif (SP_kind(you) == 5) then  ! left bottom wall
            VM_x(me,1) = 2.0d0*(WL_thick + 0.0d0) - SP_x(you,1)
            VM_x(me,2) = 2.0d0*(WL_thick + 0.0d0) - SP_x(you,2)
            VM_kind(me) = 5

        elseif(SP_kind(you) == 6) then  ! right bottom wall
            VM_x(me,1) = 2.0d0*(WL_thick + 1.0d0) - SP_x(you,1)
            VM_x(me,2) = 2.0d0*(WL_thick + 0.0d0) - SP_x(you,2)
            VM_kind(me) = 6

        elseif (SP_kind(you) == 7) then  ! left top wall
            VM_x(me,1) = 2.0d0*(WL_thick + 0.0d0) - SP_x(you,1)
            VM_x(me,2) = 2.0d0*(WL_thick + 1.0d0) - SP_x(you,2)
            VM_kind(me) = 7

        elseif (SP_kind(you) == 8) then  ! right top wall
            VM_x(me,1) = 2.0d0*(WL_thick + 1.0d0) - SP_x(you,1)
            VM_x(me,2) = 2.0d0*(WL_thick + 1.0d0) - SP_x(you,2)
            VM_kind(me) = 8

        else
            stop ' [error] @set_VN.f90'
        endif

    enddo
    
end subroutine set_VM

! END !