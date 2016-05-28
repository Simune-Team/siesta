! Handlers
! ------------------------------------
!
! Example of routine die
! In your program, depending on your own usage details,
! you might need something more sophisticated
!
subroutine die(str)
 character(len=*), intent(in) :: str

 write(*,*) trim(str)
 stop
end subroutine die
!
! Handlers for alloc
! In your program, depending on your own usage details,
! you might need something more sophisticated
!
subroutine alloc_memory_event(bytes,name)
integer, intent(in) :: bytes
character(len=*), intent(in) :: name
!! write(*,*) "alloc: allocated ", bytes, "bytes for "//trim(name)
end subroutine alloc_memory_event

subroutine alloc_error_report(name,code)
character(len=*), intent(in) :: name
integer, intent(in) :: code
!! write(*,*) "alloc error: "//trim(name)
end subroutine alloc_error_report
!
! Timer interface: called at specific events
! It will include MPI calls if gridxc_init included
! the relevant option

         subroutine gridxc_timer_start(str)
         character(len=*), intent(in)  :: str
         ! print *, "Entered: "//trim(str)
         end subroutine gridxc_timer_start

         subroutine gridxc_timer_stop(str)
         character(len=*), intent(in)  :: str
         ! print *, "Exited: "//trim(str)
         end subroutine gridxc_timer_stop
!
