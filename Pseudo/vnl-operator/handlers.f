!--------------------------------------------------
! Stand-alone routine to capture error messages from
! the alloc module
!
! This functionality could be made more general,
! and use a uniform interface for all the utility
! modules developed in-house. (Let's say, call it
! 'error_report' with severity arguments, etc).
!
! Each program using the alloc module needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta:
!
!   - The tagging of the message by Node.
!   - The use of 'unit 6' as output and '0' as error.
!
! Routines using this functionality should include
! the following
!
!   subroutine alloc_error_report(str,code)
!     character(len=*), intent(in) :: str
!     integer, intent(in)          :: code
!   end subroutine alloc_error_report
!
!------------------------------------------------------

      subroutine alloc_error_report(str,code)

      character(len=*), intent(in)  :: str
      integer, intent(in)  :: code

      external :: die
      
      if (code == 0) then
        call die(str)
      else
         write(0,*) trim(str)
         write(6,*) trim(str)
      endif

      end subroutine alloc_error_report
!--------------------------------------------------
! Stand-alone routine to capture memory events from
! the alloc module
!
! Each program using the alloc module needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta
!
!   - The use of the memory_log module for reports
!
! Routines using this functionality should include
! the following
!
!   subroutine alloc_memory_event(str,code)
!     character(len=*), intent(in) :: str
!     integer, intent(in)          :: code
!   end subroutine alloc_memory_event
!
!------------------------------------------------------

      subroutine alloc_memory_event(bytes,name)

      integer, intent(in)           :: bytes
      character(len=*), intent(in)  :: name

      ! Do nothing

      end subroutine alloc_memory_event

