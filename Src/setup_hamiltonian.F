      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( iscf )

      USE siesta_options
      use sparse_matrices, only: H_kin, H_vkb, H0
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
      use m_spin,         only: nspin, SpOrb
      use m_spin,         only: h_spin_dim, e_spin_dim, spinor_dim
      use spinorbit,      only: spinorb
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

c      stop 'Stopping in setup_hamiltonian 1....'   


      ! Initialize Hamiltonian with the non-SCF part...................

CC RC      do ispin = 1,nspin  ! CC RC Check this for osSO 
CC RC        if (ispin .le. 2) then
CC RC          do io = 1,maxnh
CC RC            H(io,ispin) = H_kin(io) + H_vkb(io)
CC RC          enddo
CC RC        else
CC RC          do io = 1,maxnh
CC RC            H(io,ispin) = 0.0_dp
CC RC         enddo
CC RC       endif
CC RC      enddo          

      ! Initialize Hamiltonian with the non-SCF part...................

      if (SpOrb) then
       H(:,:) = H0(:,:)

c       write(6,'(a)') 'kkkkkkkkkkkkkk'
c        write(6,'(a)') ' setup_hamiltonian: Writing H elem '
c        do io = 1,maxnh
c          write(6,'(a,i5,8f12.6)') 'ind/H0(io,1:8) = ',
c     .             io, H0(io,1:8)
c        enddo
c       write(6,'(a)') 'kkkkkkkkkkkkkk'
      else
       H(:,:) = 0.0_dp
       do ispin = 1,spinor_dim
        H(:,ispin) = H_kin(:,1) + H_vkb(:,1)! Thios could be H0(:,1)
       enddo
      endif

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
      do ispin = 1,min(nspin,2)  ! CC RC This is ok, =spinor_dim
        do io = 1,maxnh
c             write(6,'(a,2i5,2f12.6)') 'ind/Dscf(ind,ispin) = ',
c     .             ispin,io, Dscf(io,ispin)
          Ekin = Ekin + H_kin(io,1) * Dscf(io,ispin)
          Enl = Enl + H_vkb(io,1) * Dscf(io,ispin)
        enddo
      enddo

CC RC
#ifdef MPI
! Global reduction of Ekin, Enl and Eso 
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
#endif

c      write(6,'(a,f12.6)') 'setup_hamiltonian: Ekin = ', Ekin
c      write(6,'(a,f12.6)') 'setup_hamiltonian:  Enl = ', Enl
c      write(6,'(a,f12.6)') ' setup_hamiltonian:      Ekin+Enl = ', 
c     .                       Ekin+Enl

CC RC for osSO
c      write(6,'(a,f12.6)') ' before assign zero value to Eso = ', Eso
      Eso = 0.0_dp
      if (SpOrb) then
c       write(6,'(a)') 'kkkkkkkkkkkkkk'
       do io = 1,maxnh
c        write(6,'(a,6f12.6)') 'Dscf(io,3:8) = ', Dscf(io,3:8)
c              write(6,'(a,6f12.6)') 'H0(io,3:8) = ', H0(io,3:8)
        Eso = Eso + H0(io,3)*Dscf(io,7) + H0(io,4)*Dscf(io,8)
     .            + H0(io,7)*Dscf(io,3) + H0(io,8)*Dscf(io,4)
     .            - H0(io,5)*Dscf(io,5) - H0(io,6)*Dscf(io,6)
       enddo
c       write(6,'(a)') 'kkkkkkkkkkkkkk'
#ifdef MPI
! Global reduction of Ekin, Enl and Eso 
      call globalize_sum(Eso,buffer1)
      Eso = buffer1
#endif
c       write(6,'(a,f12.6)') ' setup_hamiltonian:           Eso = ', Eso
      endif
CC RC 
      E0 = Ekin + Enl + Eso - Eions + Ena
c      write(6,'(a,f12.6)') ' setup_hamiltonian:    Eso = ',
c     .    E0-Ekin-Enl+Eions-Ena 

CC RC
c      write(6,'(a,f12.6)') ' setup_hamiltonian:  Ekin+Enl+Eso = ', 
c     .                       Ekin+Enl+Eso 
C          stop 'setup_hamiltonian: Stopping after the second step'

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
  
c      write(6,'(a,2i3)') 'nspin/h_spin_dim = ', nspin, h_spin_dim 
c      stop 'Stopping before call to dhscf...'

      call dhscf( nspin, h_spin_dim, SpOrb, no_s, iaorb, iphorb, no_l,
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
