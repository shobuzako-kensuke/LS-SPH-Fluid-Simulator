subroutine set_xyz_kind
    !=========================!
    !  module                 ! 
    !=========================!
    use input
    use global_variables

    !=========================!
    !  local variables        ! 
    !=========================!
    implicit none
    integer :: N_x_fill, N_y_fill, N_sys_fill, i
    real(8), allocatable :: SP_xyz_sys(:,:), SP_xyz_inn(:,:), SP_xyz_out(:,:)
    
    !==============================!
    !  fill system with particles  !
    !==============================!
    N_x_fill   = N_x + 2*N_WL        ! inner + outer particles along x-axis
    N_y_fill   = N_y + 2*N_WL        ! inner + outer particles along y-axis
    N_sys_fill = N_x_fill * N_y_fill ! all particles if system packed

    allocate(SP_xyz_sys(N_sys_fill, 2)) ! (name, (x,y))
    allocate(SP_xyz_inn(N_sys_fill, 2))
    allocate(SP_xyz_out(N_sys_fill, 2))

    SP_xyz_sys = 0.0d0 ! zero clear
    SP_xyz_inn = 0.0d0
    SP_xyz_out = 0.0d0
    SP_xyz_sys(1,:) = Delta_x / 2.0d0

    do i = 2, N_sys_fill
        if (i <= N_x_fill) then
            SP_xyz_sys(i,1) = SP_xyz_sys(i-1,1) + Delta_x
            SP_xyz_sys(i,2) = SP_xyz_sys(i-1,2)
        else
            SP_xyz_sys(i,1) = SP_xyz_sys(i-N_x_fill, 1)
            SP_xyz_sys(i,2) = SP_xyz_sys(i-N_x_fill, 2) + Delta_x
        endif
    enddo

    !=========================!
    !  categorize             !
    !=========================!
    N_inn = 0 ! zero clear
    N_out = 0

    do i = 1, N_sys_fill
        !=========================!
        !  find inner particles   !
        !=========================!
        if ((WL_thick + 0.0d0  < SP_xyz_sys(i,1)) .and. &
            (WL_thick + width  > SP_xyz_sys(i,1)) .and. &
            (WL_thick + 0.0d0  < SP_xyz_sys(i,2)) .and. &
            (WL_thick + height > SP_xyz_sys(i,2))) then
            
            N_inn = N_inn + 1 ! count the number of inner particles
            SP_xyz_inn(N_inn,:) = SP_xyz_sys(i,:)

        !==============================!
        !  find wall particles         !
        !==============================!
        else
            N_out = N_out + 1 ! count the number of outer particles
            SP_xyz_out(N_out,:) = SP_xyz_sys(i,:)
        endif
    enddo

    !=========================!
    !  allocate & sort        !
    !=========================!
    N_sys = N_inn + N_out

    allocate(SP_xyz (N_sys,2)) ! (name, (x,y))
    allocate(SP_kind(N_sys))   ! (name)
    SP_xyz  = 0.0d0            ! zero clear
    SP_kind = 0                ! zero clear

    ! inner particles >> SP_kind = 0
    SP_xyz (1:N_inn, :) = SP_xyz_inn(1:N_inn, :)
    SP_kind(1:N_inn   ) = 0

    ! wall particles >> SP_kind = 10 (tmp)
    SP_xyz ((N_inn+1):, :) = SP_xyz_out(1:N_out, :)
    SP_kind((N_inn+1):   ) = 10
    
    !===========================!
    !  identify wall particles  !
    !===========================!
    do i = N_inn+1, N_sys

        ! bottom wall >> SP_kind = (1,5,6)
        if (SP_xyz(i,2) < WL_thick) then
            ! left bottom
            if (SP_xyz(i,1) < WL_thick) then
                SP_kind(i) = 5

            ! right bottom
            elseif (SP_xyz(i,1) > WL_thick + width) then
                SP_kind(i) = 6
            
            ! bottom
            else
                SP_kind(i) = 1
            endif
        
        ! top wall >> SP_kind = (2,7,8)
        elseif (SP_xyz(i,2) > WL_thick + height) then
            ! left top
            if (SP_xyz(i,1) < WL_thick) then
                SP_kind(i) = 7

            ! right top
            elseif (SP_xyz(i,1) > WL_thick + width) then
                SP_kind(i) = 8
            
            ! top
            else
                SP_kind(i) = 2
            endif

        ! left wall >> SP_kind = 3
        elseif (SP_xyz(i,1) < WL_thick) then
            SP_kind(i) = 3

        ! right wall >> SP_kind = 4
        elseif (SP_xyz(i,1) > WL_thick + width) then
            SP_kind(i) = 4

        else
            stop ' [error] SP_kind is not correct @ set_xyz_kind'
        endif
    enddo

    !=========================!
    !  deallocate             !
    !=========================!
    deallocate(SP_xyz_sys, SP_xyz_out, SP_xyz_inn)
    
end subroutine set_xyz_kind

! END !