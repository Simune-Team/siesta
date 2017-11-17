! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine atm_transfer()

      use atm_types, only: maxnorbs, nspecies
      use atm_types, only: species, species_info

      use radial
      use atmparams, only:NTBMAX
!----------------------------------------------------------------
      use old_atmfuncs, only: nsmax
!
!     old_atmfuncs arrays
!
      use old_atmfuncs, only: tabpol, table, tab2
      use old_atmfuncs, only: tab2pol
      use old_atmfuncs, only: qtb, slfe
      use old_atmfuncs, only: lmxosave, npolorbsave
      use old_atmfuncs, only: nzetasave, nsemicsave
!
!     old_atmfuncs procedures
!
      use old_atmfuncs, only: labelfis, izofis, zvalfis
      use old_atmfuncs, only: massfis, lomaxfis, nofis
      use old_atmfuncs, only: cnfigfio, lofio, mofio
      use old_atmfuncs, only: atmpopfio, rcut

!----------------------------------------------------------------
      use ldau_specs,     only: populate_species_info_ldau

      use periodic_table, only: symbol
      use sys,            only: die

      implicit none

      type(species_info), pointer        :: spp
      type(rad_func), pointer            :: op
      type(rad_func), pointer            :: pp

      integer is, io, i , n, ntot, l, m
      integer max_norbnl, nsm, izeta, j
      integer norb, indx, ipol, num_normal, num_pol

      integer, dimension(maxnorbs) :: index_normal, z_normal,
     $     nsm_normal, index_pol, z_pol, nsm_pol
      integer, dimension(maxnorbs) :: mark  ! overkill

!     Fill in main structures for new atmfuncs -----------

      nspecies = nsmax           ! From old_atmfuncs

!-----------------------------------------------------------

      max_norbnl = 0
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

!        Orbitals done in Basis_gen (for SPLIT and vanilla POLgen)
!
!        KB projectors, done in kbgen
!        chlocal, vlocal, vna, core, done in atom

      enddo

      call populate_species_info_ldau

      end subroutine atm_transfer







