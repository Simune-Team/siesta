! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      module sys
!
!     Termination and messaging routines, MPI aware
!     This module is being progressively emptied out      
!     of functionality, relying instead in properly
!     adapted external functions
      
      implicit none

      public :: die      ! Prints an error message and terminates the program
      public :: message  ! Prints a message string if node==0

      private

      CONTAINS

      subroutine message(level,str)
      
      character(len=*), intent(in)  :: level
      ! One of INFO, WARNING, FATAL
      character(len=*), intent(in)  :: str

      interface
         subroutine message_external(level,str)
         character(len=*), intent(in) :: level
         character(len=*), intent(in) :: str
         end subroutine message_external
      end interface

      call die_external(level,str)

      end subroutine message
!
!--------------------------------------------------
      subroutine die(str)
      use siesta_cml
      character(len=*), intent(in), optional   :: str
      interface
         subroutine die_external(str)
         character(len=*), intent(in) :: str
         end subroutine die_external
      end interface

      if (present(str)) then
         call die_external(str)
      else
         call die_external("sys::die: Program terminated")
      endif
      end subroutine die

      end module sys

!
!     Special versions for use in sys.F
!      
      subroutine die_external(str)
      character(len=*), intent(in)  :: str
      external :: die
      call die(str)
      end subroutine die_external

      subroutine message_external(level,str)
      character(len=*), intent(in)  :: level
      character(len=*), intent(in)  :: str
      external :: message
      call message(level,str)
      end subroutine message_external

