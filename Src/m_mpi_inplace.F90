!     
!     Copyright (C) 1996-2016	The SIESTA group
!     This file is distributed under the terms of the
!     GNU General Public License: see COPYING in the top directory
!     or http://www.gnu.org/copyleft/gpl.txt.
!     See Docs/Contributors.txt for a list of contributors.
!     
!> Provides interfaces for in-place all-reduce operations
!> This operation is not supported by the old-style interfaces
!> in the mpi_siesta module. Hence, the generic 'mpi' module
!> is used instead, and the communicator argument is mandatory,
!> as we cannot access the (possibly redefined) mpi_comm_world
!> as a default.
!> Note that (when implemented) a scalar-integer interface in this module
!> will not be distinguishable from that in m_mpi_utils.

module m_mpi_inplace
      use precision, only: dp, sp
#ifdef MPI
      use mpi
      implicit none
      integer, private :: MPIerror
#else
      implicit none
#endif
      public :: globalize_sum
      private
      
      interface globalize_sum
       module procedure globalize_sum_inplace_dp
       module procedure globalize_sum_inplace_v_dp
       module procedure globalize_sum_inplace_vv_dp
      end interface

      CONTAINS

!--------------------------------------------------------------
!     In-place versions of 'dp' interfaces
      
      subroutine Globalize_sum_inplace_dp(var,comm)
      real(dp), intent(inout) :: var
      integer, intent(in)     :: comm

#ifdef MPI
      call MPI_AllReduce(MPI_IN_PLACE,var,1,MPI_double_precision, &
                         MPI_sum,comm,MPIerror)
#else
      ! do nothing
#endif
      end subroutine Globalize_sum_inplace_dp

      subroutine Globalize_sum_inplace_v_dp(var,comm)
      real(dp), intent(inout), dimension(:)  :: var
      integer, intent(in)                    :: comm

#ifdef MPI
      call MPI_AllReduce(MPI_IN_PLACE,var(1),size(var),MPI_double_precision, &
                         MPI_sum,comm,MPIerror)
#else
      ! do nothing
#endif
      end subroutine Globalize_sum_inplace_v_dp

      subroutine Globalize_sum_inplace_vv_dp(var,comm)
      real(dp), intent(inout), dimension(:,:)  :: var
      integer, intent(in)                      :: comm

#ifdef MPI
      call MPI_AllReduce(MPI_IN_PLACE,var(1,1),size(var),MPI_double_precision, &
                        MPI_sum,comm,MPIerror)
#else
      ! do nothing
#endif
      end subroutine Globalize_sum_inplace_vv_dp

    end module m_mpi_inplace

