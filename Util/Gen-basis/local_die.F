! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!--------------------------------------------------
! Stand-alone 'die' routine for use by libraries and
! low-level modules.
!
! Each program using the module or library needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta:
!
!   - The need to have the pxf functionality.
!   - The use of 'unit 6' as output.
!
! Routines using this functionality should include
! the following
!
!     interface
!      subroutine die(str)
!      character(len=*), intent(in)  :: str
!      end subroutine die
!     end interface
!
!------------------------------------------------------

      subroutine die(str)

      character(len=*), intent(in)  :: str

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)
         call pxfflush(6)
         call pxfflush(0)
      call pxfabort()
      end subroutine die

      subroutine psml_die(str)
      character(len=*), intent(in)  :: str
      call die(str)
      end subroutine psml_die
