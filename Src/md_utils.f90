! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module md_utils

  use precision, only:  dp

  integer, save, public  :: md_io_unit = 99
  logical, save, public  :: restart_with_mdinfo  = .false.
  logical, save,  :: file_opened  = .false.

  CONTAINS

    subroutine add_to_md_file(natoms,xa,va,cell,vcell,nose,nosedot)

      integer, intent(in)                                  :: natoms
      real(dp), dimension(3,natoms), intent(in)            :: xa, va
      real(dp), dimension(3,3), intent(in), optional       :: cell, vcell
      real(dp), intent(in), optional                       :: nose, nosedot

!     Here we can save x, xa, va for MD  (experimental)
!
        write(99,*) natoms
        do ia = 1,natoms
          write(99,'(i4,3f14.8,3x,3f14.8)')
     .      iza(ia),(xa(i,ia),i=1,3),(va(i,ia),i=1,3)
        enddo

        have_nose = present(nose)
        have_nosedot = present(nosedot)

        if (have_nose) then
           write(99,*) nose
        endif

        end subroutine add_to_md_file
