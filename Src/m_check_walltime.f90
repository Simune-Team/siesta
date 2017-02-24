module m_check_walltime

  
  integer, public :: WALLTIME_MAX
  integer, public :: WALLTIME_WARNING

  integer, private, parameter :: dp = selected_real_kind(10,100)
  
  contains
    subroutine check_walltime(time_is_up)
      use m_walltime, only: wall_time

      logical, intent(out) :: time_is_up
      real(dp) :: t

      call wall_time(t)
      time_is_up = (t > WALLTIME_WARNING)
    end subroutine check_walltime

  end module m_check_walltime
  
