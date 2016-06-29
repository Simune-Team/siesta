! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_normalize_dm
      private
      public :: normalize_dm

      CONTAINS

      subroutine normalize_dm( first )
      use precision, only: dp
      use sparse_matrices, only: Dscf, Escf, maxnh, S
      use siesta_options,  only: dm_normalization_tol
      use siesta_options,  only: normalize_dm_during_scf
      use atomlist, only: qtot
      use parallel, only: IOnode
      use sys,      only: die
      use m_spin,   only: nspin
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
      implicit none
      logical :: first
      integer :: io, ispin
      real(dp):: qsol  ! Total unnormalized electron charge
      character(len=132) :: msg
#ifdef MPI
      real(dp):: buffer1
#endif

!     Normalize density matrix to exact charge
      qsol = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,maxnh
          qsol = qsol + Dscf(io,ispin) * s(io)
        enddo
      enddo
#ifdef MPI
      call globalize_sum(qsol,buffer1)
      qsol = buffer1
#endif

      ! We might consider using abs(qsol-qtot) instead of the fractional accuracy.
      ! For large numbers of electrons one might have a significant net charge
      ! before the fractional tolerance is reached.

      ! Two degrees of reporting:
      if (first) then
         ! If 1st scf step, it might be that a (reused) DM was not
         ! normalized properly
         if (abs(qsol/qtot-1._dp) > 1.0e-3_dp) then
            if (IOnode) then
               write(6,'(a,2f20.8)')
     &              'Note: For starting DM, Qtot, Tr[D*S] =', qtot, qsol
            endif
         endif
      else
         ! In later steps, the lack of normalization is more serious
         ! The tolerance is tighter by default
         if (abs(qsol/qtot-1._dp) > dm_normalization_tol) then
            write(msg,'(a,2f20.8)')
     $           'Bad DM normalization: Qtot, Tr[D*S] =', qtot, qsol
            call die(msg)
         endif
      endif

      ! Actually normalize only if not disallowed by the user
      !
      if (normalize_dm_during_scf) then
         call normalize()
      else
         if (IOnode) then
            write(6,'(a,g15.6)')
     &           'Note: Not normalizing DM: log(Tr[D*S]/Qtot -1) =',
     &           log10(abs((qsol/qtot)-1.0_dp))
         endif
      endif

      CONTAINS

      subroutine normalize()
      Dscf(1:maxnh,1:nspin) = Dscf(1:maxnh,1:nspin) * qtot/qsol
      Escf(1:maxnh,1:nspin) = Escf(1:maxnh,1:nspin) * qtot/qsol
      end subroutine normalize
      

      END subroutine normalize_dm

      END module m_normalize_dm




