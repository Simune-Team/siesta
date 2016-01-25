! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_compute_max_diff
      public :: compute_max_diff
      CONTAINS
      subroutine compute_max_diff(Xin,Xout,MaxDiff)

      use precision,         only: dp
      use sparse_matrices,   only: numh, listhptr
      use m_spin,            only: nspin
      use atomlist,          only: no_l
#ifdef MPI
      use m_mpi_utils,           only: globalize_max
#endif

      real(dp), intent(in)    :: Xin(:,:), Xout(:,:)
      real(dp), intent(out)   :: MaxDiff

      integer   :: is, i, ind, in
#ifdef MPI
      real(dp)              :: buffer1
#endif

      MaxDiff = 0.0_dp
      do is = 1,nspin
         do i = 1,no_l
            do in = 1,numh(i)
               ind = listhptr(i) + in
               MaxDiff = max(MaxDiff, abs(Xout(ind,is) - Xin(ind,is)))
            enddo
         enddo
      enddo

#ifdef MPI
!     Ensure that MaxDiff is the same on all nodes 
      call globalize_max(MaxDiff,buffer1)
      MaxDiff = buffer1
#endif
      end subroutine compute_max_diff
      end module m_compute_max_diff
