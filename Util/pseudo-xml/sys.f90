! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module sys
!
!     Termination and messaging routines, MPI aware
!
      public :: die, bye, message

      CONTAINS

      subroutine message(str)

      character(len=*), intent(in), optional   :: str


         if (present(str)) then
            write(6,'(a)') trim(str)
         endif

      end subroutine message
!
!--------------------------------------------------
      subroutine die(str)

      character(len=*), intent(in), optional   :: str

         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
      STOP
      end subroutine die

!---------------------------------------------------------
      subroutine bye(str)

      character(len=*), intent(in), optional   :: str

         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Requested End of Run. Bye!!'
!!                                       endif

      stop

      end subroutine bye

      end module sys

