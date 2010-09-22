      module m_setup_H0
      private
      public :: setup_H0
      CONTAINS

      subroutine setup_H0()
      
C     Computes non-self-consistent part of the Hamiltonian
C     and initializes data structures on the grid.
      
      USE m_normalize_dm, only: normalize_dm
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     &                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     &                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use m_dhscf,      only: dhscf_init
      use m_stress
      use parallel,     only: Node
      use m_energies
      use m_steps, only: istp
      use m_ntm
      use m_spin,         only: nspin
      use m_dipol
      use m_forces,          only: fa
      use alloc, only: re_alloc, de_alloc
      use m_gamma
      use sys, only: die
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
      use mpi_siesta
#endif
#ifdef TRANSIESTA
!     TSS Begin
      use m_ts_options, only: TSrun
!     TSS End
#endif /* TRANSIESTA */

      implicit none
      real(dp) :: stressl(3,3), dummy(1,1)
      integer  :: ia, io, is
      real(dp) :: g2max
#ifdef MPI
      real(dp) :: buffer1
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE setup_H0' )
#endif
!     Normalize density matrix to exact charge
#ifdef TRANSIESTA      
      if (.not.TSrun ) call normalize_dm( first=.true. )
#else
      call normalize_dm( first=.true. )
#endif  /* TRANSIESTA */

C     Initialize Hamiltonian
      H = 0.0_dp

C     Self-energy of isolated ions
      Eions = 0.0_dp
      do ia = 1,na_u
        is = isa(ia)
        Eions = Eions + uion(is)
      enddo

C     In these routines, add a flag to tell them NOT to compute
C     forces and stresses in this first pass.   ********

C     Neutral-atom: energy, forces and stress ............................
C     First time for energy, last time for forces
      call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &           isa, Ena, fa, stress, forces_and_stress=.false.)
      call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, DEna, fa, stress, forces_and_stress=.false.) 
      Ena = Ena + DEna
      
C     Kinetic: energy, forces, stress and matrix elements .................
      call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &             maxnh, maxnh, lasto, iphorb, isa, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             min(nspin,2), Dscf, Ekin, dummy, stressl, H,
     &             forces_and_stress=.false.) 
#ifdef MPI
C     Global reduction of energy terms
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
#endif

C     Non-local-pseudop: energy, forces, stress and matrix elements .......
      call nlefsm(scell, na_u, na_s, isa, xa, indxua, 
     &            maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     &            numh, listhptr, listh, numh, listhptr, listh, 
     &            min(nspin,2), Dscf, Enl, dummy, stressl, H,
     &            forces_and_stress=.false.)

#ifdef MPI
C     Global reduction of energy terms
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
#endif

C     Save partial Hamiltonian (non-SCF part) ......................
      do io = 1,nh
        H0(io) = H(io,1)
      enddo

C     Non-SCF part of total energy .......................................
      E0 = -Eions + Ena + Ekin + Enl

C     Allocate all the structures needed for grid work (formerly
C     in the first section of dhscf  *******

C     Output memory use after setting up the structures
      if (istp.eq.1) call printmemory( 6, 0 )

C     Metadynamics energy
      if (lMetaForce) then
        call meta(xa,na_u,ucell,Emeta,fa,stress,.false.,.false.)
      endif

C     Add on force field contribution to energy
      call twobody( na_u,xa,isa,ucell,Emm,
     &              ifa=0,fa=fa,istr=0,stress=stress)

C     This will take care of possible changes to the mesh and atomic-related
C     mesh structures for geometry changes
      g2max = g2cut
      call dhscf_init( nspin, no_s, iaorb, iphorb,
     &                 no_l, no_u, na_u, na_s,
     &                 isa, xa, indxua, ucell,
     &                 mscell, G2max, ntm,
     &                 maxnh, numh, listhptr, listh, datm,
     &                 dummy, stressl)

#ifdef DEBUG
      call write_debug( '  POS setup_H0' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine setup_H0
      END module m_setup_H0
