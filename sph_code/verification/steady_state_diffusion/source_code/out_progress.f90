subroutine out_progress
    !=========================!
    !  module                 !
    !=========================!
    use input
    use global_variables

    !=========================!
    !  local variables        !
    !=========================!
    implicit none
    character(len=1  ) :: esc = achar(27)
    character(len=999) :: progress, residual, progress_time

    !=========================!
    !  progress bar           !
    !=========================!
    write(progress, *) step
    write(residual, '(es12.3)') residual_max
    write(progress_time, '(f12.3)') tmp_time - start_time

    write(*,'(a)', advance='no') esc//'[1A'//esc//'[2K'
    write(*,'(a)', advance='no') esc//'[1A'//esc//'[2K'
    write(*,'(a)', advance='no') esc//'[1A'//esc//'[2K'

    write(*,'(1x,a)') '[message] step     : '//trim(adjustl(progress))
    write(*,'(1x,a)') '[message] residual : '//trim(adjustl(residual))
    write(*,'(1x,a)') '[message] CPU time : '//trim(adjustl(progress_time))//' [s]'
    
end subroutine out_progress

! END !
