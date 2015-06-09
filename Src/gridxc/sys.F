!!@LICENSE

      module sys

      ! This module contains interfaces to routines
      ! that must be provided as external by the user
      ! of the library

      interface
         subroutine die(str)
         character(len=*), intent(in)  :: str
         end subroutine die
      end interface
      interface
         subroutine gridxc_timer_start(str)
         character(len=*), intent(in)  :: str
         end subroutine gridxc_timer_start
      end interface
      interface
         subroutine gridxc_timer_stop(str)
         character(len=*), intent(in)  :: str
         end subroutine gridxc_timer_stop
      end interface
      interface
         subroutine memory_accouting_event(name,bytes)
         character(len=*), intent(in)  :: name
         integer, intent(in)           :: bytes
         end subroutine memory_accouting_event
      end interface
      end module sys

