      MODULE m_normalize_dm
      private
      public :: normalize_dm

      CONTAINS

      subroutine normalize_dm( first )
      USE SIESTA_TODO
      use precision
      use sparse_matrices
      use atomlist, only: qtot
      use parallel, only: IOnode
      use m_spin,   only: nspin

      implicit none

      logical :: first
      integer :: io, ispin
      real(dp):: qsol       ! Total unnormalized electron charge

      ! Normalize density matrix to exact charge
      qsol = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,nh
          qsol = qsol + Dscf(io,ispin) * s(io)
        enddo
      enddo
#     ifdef MPI
        call globalize_sum(qsol,buffer1)
        qsol = buffer1
#     endif
      if (IOnode) then
        if (.not.first .and.
     .     abs(qsol/qtot-1._dp).gt.1.d-2) write(6,'(a,2f15.6)')
     .    'siesta: WARNING: Qtot, Tr[D*S] =', qtot, qsol
      endif
      do ispin = 1,nspin
        do io = 1,nh
          Dscf(io,ispin) = Dscf(io,ispin) * qtot/qsol
          Escf(io,ispin) = Escf(io,ispin) * qtot/qsol
        enddo
      enddo

      END subroutine normalize_dm

      END module m_normalize_dm




