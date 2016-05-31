! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_mixer
      private
      public :: mixer

      CONTAINS
 
      subroutine mixer( iscf )

      use precision,        only: dp
      use siesta_options,   only: mixH, mix_scf_first, mix_spinor
      use siesta_options,   only: mullipop, muldeb
      use sparse_matrices,  only: Dold, Dscf, Hold, H, S
      use sparse_matrices,  only: maxnh, numh, listhptr, listh
      use siesta_geom,      only: na_u, isa
      use atomlist,         only: iaorb, iphorb, lasto, no_u, no_l

      use m_mixing_scf, only: scf_mix
      use m_mixing, only: mixing
      
      use m_spin,           only: nspin
      use parallel,         only: IONode

      implicit none

      real(dp), pointer :: Xin(:,:), Xout(:,:)
      real(dp), allocatable :: F(:,:)

      integer,  intent(in)  :: iscf

      real(dp)              :: dummy_dmax   ! Now local variable
      integer               :: iiscf
      logical               :: mmix    ! Same as mix. Used in pulayx

      external :: mulliken

!-------------------------------------------------------------------- BEGIN

      call timer( 'MIXER', 1 )

      if ( mixH ) then
         ! Mix Hamiltonian
         Xin => Hold
         Xout => H
      else
         ! Mix density matrix
         Xin => Dold
         Xout => Dscf
      end if


      ! Create residual function to minimize
      allocate(F(maxnh,nspin))
!$OMP parallel workshare default(shared)
      F = Xout - Xin
!$OMP end parallel workshare

      ! Call mixing routine
      ! Xin contains the input element, F contains
      !  F = Xout - Xin
      ! Upon exit Xout contains the mixed quantity
      if ( mix_spinor ) then
         call mixing( scf_mix, iscf, maxnh, nspin, Xin, F, Xout,
     &        ncoeff=min(2,nspin))
      else
         call mixing( scf_mix, iscf, maxnh, nspin, Xin, F, Xout)
      end if
      
      ! Correctly handle mixed quantity
      ! move over mixed quantity to "input"
      if ( mix_scf_first ) then
         
         ! always allow mixing
         Xin = Xout
         
      else if ( iscf == 1 ) then
         
         ! We are not allowed to mix the first SCF
         ! Remember that Xout contains the MIXED
         ! quantity and F is the current residual
         ! hence:
         !   Xout = Xin + F
!$OMP parallel workshare default(shared)
         Xin = Xin + F
         ! we require that Xin == Xout, so this is
         ! required:
         Xout = Xin
!$OMP end parallel workshare

      end if

      deallocate(F)

      ! Print populations at each SCF step, if requested
      ! Note that this is after mixing, which is not
      ! entirely correct. It should be moved to the top,
      ! or done somewhere else.

      if (muldeb) then 
         if (IONode)
     &        write (6,"(/a)")
     &        'siesta: Mulliken populations after mixing'
        call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     &                 numh, listhptr, listh, S, Dscf, isa,
     &                 lasto, iaorb, iphorb )
      endif

      call timer( 'MIXER', 2 )

!-------------------------------------------------------- END
      END subroutine mixer

      End MODULE m_mixer


