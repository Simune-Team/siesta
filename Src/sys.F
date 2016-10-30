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

      public :: die      ! Prints an error message and calls MPI_Abort
      public :: message  ! Prints a message string if node==0

      private

      CONTAINS

      subroutine message(str)

#ifdef MPI
      use mpi_siesta
#endif
      
      character(len=*), intent(in), optional   :: str

      integer :: Node
#ifdef MPI
      integer :: MPIerror
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
      endif

      end subroutine message
!
!--------------------------------------------------
      subroutine die(str)
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
!     Special version for use in sys.F
!      
      subroutine die_external(str)
      character(len=*), intent(in)  :: str
      external :: die
      call die(str)
      end subroutine die_external

