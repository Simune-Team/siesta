! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!--------------------------------------------------


      subroutine timer(str,iopt)
      character(len=*), intent(in)  :: str
      integer, intent(in)  :: iopt
      ! do nothing
      end subroutine timer

      
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
      use memory_log, only: memory_event

      integer, intent(in)           :: bytes
      character(len=*), intent(in)  :: name

      call memory_event(bytes,name)

      end subroutine alloc_memory_event

! The PSML library calls a "die" routine when it encounters an
! error. This routine should take care of carrying out any needed
! cleaning and terminating the program.  As the details would vary with
! the client program, each program has to provide its own.
! 

      subroutine psml_die(str)
      character(len=*), intent(in) :: str

      write(0,"(a)") str
      STOP

      end subroutine psml_die
