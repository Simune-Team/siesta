! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( iscf )

      USE siesta_options
      use sparse_matrices, only: H_kin, H_vkb
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H, S, Hold
      use sparse_matrices, only: Dscf, Escf, xijo
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_dhscf,      only: dhscf
      use m_stress
      use m_energies
      use parallel, only: Node
      use m_steps, only: istp
      use m_ntm
      use m_spin,         only: nspin
      use m_dipol
      use alloc, only: re_alloc, de_alloc
      use m_gamma
      use m_hsx, only: write_hsx
      use sys, only: die, bye
      use m_partial_charges, only: want_partial_charges
      use files, only : filesOut_t    ! derived type for output file names
      use m_rhog, only: rhog_in, rhog
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
#ifdef TRANSIESTA
      use m_ts_global_vars, only: TSrun
#endif /* TRANSIESTA */
      implicit none
      integer, intent(in) :: iscf
      real(dp)            :: stressl(3,3)
      real(dp), pointer   :: fal(:,:)   ! Local-node part of atomic F
#ifdef MPI
      real(dp)            :: buffer1
#endif
      integer             :: io, is, ispin
      integer             :: ifa     ! Calc. forces?      0=>no, 1=>yes
      integer             :: istr    ! Calc. stress?      0=>no, 1=>yes
      integer             :: ihmat   ! Calc. hamiltonian? 0=>no, 1=>yes
      real(dp)            :: g2max
      type(filesOut_t)    :: filesOut  ! blank output file names
      logical             :: use_rhog_in

!------------------------------------------------------------------------- BEGIN

      call timer('setup_H',1)

! Save present H matrix
      if(mixH) then
         Hold(1:maxnh,1:nspin) = H(1:maxnh,1:nspin)
      end if


      ! Initialize Hamiltonian with the non-SCF part...................

      do ispin = 1,nspin
        if (ispin .le. 2) then
          do io = 1,maxnh
            H(io,ispin) = H_kin(io) + H_vkb(io)
          enddo
        else
          do io = 1,maxnh
            H(io,ispin) = 0.0_dp
          enddo
        endif
      enddo          
! ..................

! Non-SCF part of total energy .......................................
! Note that these will be "impure" for a mixed Dscf

! If mixing the charge, Dscf is the previous step's DM_out. Since
! the "scf" components of the energy are computed with the (mixed)
! charge, this introduces an inconsistency. In this case the energies
! coming out of this routine need to be corrected.
! 
      Ekin = 0.0_dp
      Enl = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,maxnh
          Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
          Enl = Enl + H_vkb(io) * Dscf(io,ispin)
        enddo
      enddo
#ifdef MPI
! Global reduction of Ekin and Enl
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
#endif
      E0 = Ekin + Enl - Eions + Ena


! Add SCF contribution to energy and matrix elements ..................
      g2max = g2cut

      nullify(fal)
      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian' )

      ifa  = 0
      istr = 0
      ihmat = 1
      if ((hirshpop .or. voropop)
     $     .and. partial_charges_at_every_scf_step) then
         want_partial_charges = .true.
      endif
      use_rhog_in =  (mix_charge .and. iscf > 1)
      call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua, 
     .            ntm, ifa, istr, ihmat, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, stress, fal, stressl,
     .            use_rhog_in)

      ! This statement will apply to iscf = 1, for example, when
      ! we do not use rhog_in. Rhog here is always the charge used to
      ! build H, that is, rhog_in.
      if (mix_charge) rhog_in = rhog

      want_partial_charges = .false.
      call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 

!  It is wasteful to write over and over H and S, as there are
!  no different files.
! Save Hamiltonian and overlap matrices ............................
! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file

        if (savehs .or. write_coop) then
          call write_hsx( gamma, no_u, no_s, nspin, indxuo,
     &                    maxnh, numh, listhptr, listh, H, S, qtot,
     &                    temp, xijo)
        endif

      call timer('setup_H',2)

      if (h_setup_only) then
         call timer( 'all', 2 )  ! New call to close the tree
         call timer( 'all', 3 )
         call bye("H-Setup-Only requested")
      endif

!------------------------------------------------------------------------- END
      END subroutine setup_hamiltonian
      END module m_setup_hamiltonian
