! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine atm_transfer()
!----------------------------------------------------------------
      use ldau_specs,     only: populate_species_info_ldau

      implicit none
!-----------------------------------------------------------

      ! All work done in atom now

      call populate_species_info_ldau

      end subroutine atm_transfer







