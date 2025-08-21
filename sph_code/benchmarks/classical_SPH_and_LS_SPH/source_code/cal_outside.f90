subroutine cal_outside
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
    integer :: me

    !=========================!
    !  if out of range        ! 
    !=========================!
    !$ call omp_set_dynamic(.false.)
    !$ call omp_set_num_threads(N_threads)
    !$OMP parallel default(none) &
    !$OMP private(me) &
    !$OMP shared(N_inn, WL_thick, SP_xyz, SP_uvw)
    !$OMP do
    do me = 1, N_inn
        !=========================!
        !  left                   ! 
        !=========================!
        if (SP_xyz(me,1) < WL_thick) then            
            SP_xyz(me,1) = WL_thick
            SP_uvw(me,1) = 0.0d0
        endif

        !=========================!
        !  right                  ! 
        !=========================!
        if (SP_xyz(me,1) > WL_thick + width) then
            SP_xyz(me,1) = WL_thick + width
            SP_uvw(me,1) = 0.0d0
        endif

        !=========================!
        !  bottom                 ! 
        !=========================!
        if (SP_xyz(me,2) < WL_thick) then
            SP_xyz(me,2) = WL_thick
            SP_uvw(me,2) = 0.0d0
        endif

        !=========================!
        !  top                    ! 
        !=========================!
        if (SP_xyz(me,2) > WL_thick + height) then
            SP_xyz(me,2) = WL_thick + height
            SP_uvw(me,2) = 0.0d0
        endif

    enddo
    !$OMP enddo
    !$OMP end parallel

end subroutine cal_outside

! END !