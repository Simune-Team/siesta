! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine atm_transfer()

      use atm_types, only: nspecies
      use atm_types, only: species, species_info

!----------------------------------------------------------------
      use old_atmfuncs, only: nsmax
!
!     old_atmfuncs arrays
!
      use old_atmfuncs, only: qtb, slfe
!
!     old_atmfuncs procedures
!
      use old_atmfuncs, only: labelfis, izofis, zvalfis
      use old_atmfuncs, only: massfis
!----------------------------------------------------------------
      use ldau_specs,     only: populate_species_info_ldau

      use periodic_table, only: symbol
      use sys,            only: die

      implicit none

      type(species_info), pointer        :: spp

      integer is

!     Fill in main structures for new atmfuncs -----------

      nspecies = nsmax           ! From old_atmfuncs
!-----------------------------------------------------------

      do is=1,nspecies
         spp => species(is)

         spp%read_from_file = .false.

         spp%label = labelfis(is)
         spp%z     = izofis(is)
         if (spp%z.eq.-100) then
            spp%symbol = 'BS'
         else
            ! The function 'symbol' knows how to deal
            ! with (ghost) synthetics
            spp%symbol = symbol(spp%z)
         endif
         spp%zval  = zvalfis(is)
         spp%mass  = massfis(is)
         spp%self_energy  = slfe(is)

!        Orbitals done directly in Basis_gen
!        KB projectors, done in kbgen
!        chlocal, vlocal, vna, core, done in atom

      enddo

      call populate_species_info_ldau

      end subroutine atm_transfer







