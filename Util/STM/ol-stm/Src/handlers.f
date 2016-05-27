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

      subroutine die(str)

      character(len=*), intent(in)  :: str

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)

         call pxfflush(6)
         call pxfflush(0)
      call pxfabort()
      end subroutine die

      subroutine timer(str,iopt)
      character(len=*), intent(in)  :: str
      integer, intent(in)  :: iopt
      ! do nothing
      end subroutine timer

      
