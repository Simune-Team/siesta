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

      use parallel, only : Node

      character(len=*), intent(in)  :: str
      integer, intent(in)  :: code

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

      write(6,'(a,i0,1x,a)') 'Node: ', Node, trim(str)
      write(0,'(a,i0,1x,a)') 'Node: ', Node, trim(str)

      ! The convention for this reporter, which can 
      ! send multi-call messages, is that code=0
      ! signals the end of the series
      ! We choose to kill the program
      !
      ! This is actually too stringent, as a given failed
      ! allocation could be handled more gracefully at
      ! the point at which it is made. But in absence of
      ! a 'stat' argument for the re_alloc and de_alloc
      ! routines, this is the only thing we can do....
      !
      if (code == 0) call die(str)

      end subroutine alloc_error_report
