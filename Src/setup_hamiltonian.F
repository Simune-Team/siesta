      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( iscf )
      USE m_normalize_dm, only: normalize_dm
      USE siesta_options
      use sparse_matrices, only: H_kin, H_vkb
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H, S
      use sparse_matrices, only: Dold, Dscf, Eold, Escf, xijo
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
      use sys, only: die
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

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE setup_hamiltonian' )
#endif
      
!     Normalize density matrix to exact charge
!     Setup_hamiltonian is called before the density is printed out...
#ifdef TRANSIESTA      
      if (.not.TSrun ) call normalize_dm( first=.false. )
#else
      ! Normalize density matrix to exact charge, except if
      ! we are just using the DM to generate the initial rho and exit
      !
      if (.not. save_initial_charge_density) then
         call normalize_dm(first=.false.)
      endif
#endif  /* TRANSIESTA */

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

      call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua, 
     .            ntm, ifa, istr, ihmat, ' ', ' ', ' ', ' ', ' ', ' ',
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, stress, fal, stressl)

      call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 


! Find entropy ........................................................
! Entropy is returned from the call to diagon. To add to the energy
! the entropy computed from the input charge, here it is assigned to the one
! of the former SCF step

!** Can we trust the Entropy if iscf == 1 ??
      Entropy = 0.0_dp
      if (isolve .eq. SOLVE_DIAGON) then
        if (istp.gt.1 .or. iscf.gt.1) then
          Entropy = Entrop
        endif
      endif

! Save present density matrix ........................................
      do is = 1,nspin
        do io = 1,maxnh
          Dold(io,is) = Dscf(io,is)
          Eold(io,is) = Escf(io,is)
        enddo
      enddo

!  It is wasteful to write over and over H and S, as there are
!  no different files.
! Save Hamiltonian and overlap matrices ............................
! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file

        if (savehs .or. write_coop) then
          call write_hsx( gamma, no_u, no_s, nspin, indxuo,
     &                    maxnh, numh, listhptr, listh, H, S, qtot,
     &                    temp, xijo)
        endif

#ifdef DEBUG
      call write_debug( '  POS setup_hamiltonian' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine setup_hamiltonian
      END module m_setup_hamiltonian
