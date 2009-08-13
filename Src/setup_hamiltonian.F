      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( iscf )
      USE m_normalize_dm, only: normalize_dm
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_dhscf,      only: dhscf
      use siesta_cml
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
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 13 )
#endif
      
      ! Normalize density matrix to exact charge
      ! In the first iteration we have the D already normalized
      ! by setup_H0 ... ***
      call normalize_dm(first=.false.)

      ! Initialize Hamiltonian with the non-SCF part...................

      do ispin = 1,nspin
        if (ispin .le. 2) then
          do io = 1,nh
            H(io,ispin) = H0(io)
          enddo
        else
          do io = 1,nh
            H(io,ispin) = 0.0_dp
          enddo
        endif
      enddo          
! ..................

! Non-SCF part of total energy .......................................

      E0 = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,nh
          E0 = E0 + H0(io) * Dscf(io,ispin)
        enddo
      enddo
#ifdef MPI
! Global reduction of E0
      call globalize_sum(E0,buffer1)
      E0 = buffer1
#endif
      E0 = E0 - Eions + Ena
! ..................
!         print *, "hi antes de dhscf", Node
!         call pxfflush(6)

! Non-local-pseudop: energy, forces, stress and matrix elements .......
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

! ..................
!         print *, "hi despues de dhscf", Node
!         call pxfflush(6)

! Find entropy ........................................................
! Entropy is returned from the call to diagon. To add to the energy
! the entropy computed from the input charge, here it is assigned to the one
! of the former SCF step

!** Can we trust the Entropy if iscf == 1 ??
      Entropy = 0.0_dp
      if (isolve .eq. 0) then
        if (istp.gt.1 .or. iscf.gt.1) then
          Entropy = Entrop
        endif
      endif

! Save present density matrix ........................................
      do is = 1,nspin
        do io = 1,nh
          Dold(io,is) = Dscf(io,is)
          Eold(io,is) = Escf(io,is)
        enddo
      enddo

!  It is wasteful to write over and over H and S, as there are
!  no different files.
! Save Hamiltonian and overlap matrices ............................
      if (savehs) then
        call iohs( 'write', gamma, no_u, no_s, nspin, indxuo,
     $             maxnh, numh, listhptr, listh, H, S, qtot, temp,
     $             xijo )
      endif
        if (savehs .or. write_coop) then
          call write_hsx( gamma, no_u, no_s, nspin, indxuo,
     $               maxnh, numh, listhptr, listh, H, S, qtot, temp,
     $               xijo)
        endif


#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 0 )
#endif
!--------------------------------------------------------------------------- END
      END subroutine setup_hamiltonian
      END module m_setup_hamiltonian
