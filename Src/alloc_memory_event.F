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
      use memory_log, only: memory_event

      integer, intent(in)           :: bytes
      character(len=*), intent(in)  :: name

      call memory_event(bytes,name)

      end subroutine alloc_memory_event
