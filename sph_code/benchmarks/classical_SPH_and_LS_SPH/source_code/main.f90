!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!  program structure                                                        !
!---------------------------------------------------------------------------!
!  #1. check inputs and options >> save copy files                          !
!  #2. set inner and outer(wall) particles                                  !
!  #3. set virtual markers (VM)                                             !
!  #4. set initial state                                                    !
!  #5. set background cell                                                  !
!  #6. calculate cell number >> sort                                        !
!  #7. calculate VM >> wall                                                 !
!  #8. output initial state                                                 !
!      << time loop start >>                                                !
!           #9.  integrate EOC, EOM, EOE using 4th-order Runge-Kutta method !
!           #10. particle shifting                                          !
!           #11. output data                                                !
!      << time loop end >>                                                  !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

program main
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
    character(len=999) :: save_step
    real(8) :: omp_get_wtime

    write(*,*) '+ ------------------------------------------------------------------------ +'
    write(*,*) '[message] calculation has started'
    !$ start_time = omp_get_wtime()
    
    !=========================!
    !  check inputs           ! 
    !=========================!
    call check_compile_options          ! #1. check compile options
    call out_input                      ! #1. copy input and makefile
    call set_system_info                ! #1. cal Delta_x, h, time step

    !==========================!
    !  set particles (new)     ! 
    !==========================!
    if (trim(adjustl(read_name)) == 'new') then
        call set_xyz_kind               ! #2. set particles and sort
        call set_VM_xyz_kind            ! #3. set virtual markers
        call set_physical_quantity      ! #4. set initial state
    !=========================!
    !  set particles (recal)  ! 
    !=========================!
    else
        call set_recalculation          ! #2 & #4 read particles and initial settings
        call set_VM_xyz_kind            ! #3. set virtual markers
    endif

    !=========================!
    !  set cells              ! 
    !=========================!
    call set_background_cell            ! #5. set background cell and output
    call cal_cell_info                  ! #6. make cell info
    call cal_VM_to_WL                   ! #7. cal VM >> wall

    !=========================!
    !  output                 ! 
    !=========================!
    call out_system_info                ! #8. output variables
    if (trim(adjustl(read_name)) == 'new') then
        call out_data('0')              ! #8. output
    endif

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !  time loop start                                                          !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    do step = start_step, end_step
        call cal_Runge_Kutta            ! #9. renew xy, uv, pre, rho

#if defined(PS_1ST) || defined(PS_2ND)
        call cal_PS                     ! #10. particle shifting
        call cal_VM_to_WL               ! #10. cal VM >> wall
#endif

        !=========================!
        !  output results         ! 
        !=========================!
        if (mod(step, write_step) == 0) then
            write(save_step,*) step                  ! step >> str(step)
            call out_data(trim(adjustl(save_step)))  ! #11. output
        endif

        !=========================!
        !  progress bar           ! 
        !=========================!
        !$ tmp_time = omp_get_wtime()
        call out_progress
    enddo
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !  time loop end                                                            !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    call deallocate_array
    write(*,*) '+ ------------------------------------------------------------------------ +'
    write(*,*) '[message] calculation has finished'
    write(*,*) ''
    
end program main

! END !