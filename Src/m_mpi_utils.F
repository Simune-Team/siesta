      module m_mpi_utils

      use precision, only: dp, sp
      use sys, only: die
#ifdef MPI
      use mpi_siesta, only : MPI_Comm_World, MPI_integer, MPI_logical
      use mpi_siesta, only : MPI_Sum, MPI_Max, MPI_character
#ifdef NODAT
      use mpi_siesta, only: mpi_double_precision
#else
      use mpi_siesta, only: DAT_double
#endif
      integer, private :: MPIerror
#endif

      public :: globalize_max, globalize_sum, broadcast
      private
      interface globalize_max
        module procedure globalize_sum_dp, globalize_sum_int
      end interface
      interface globalize_sum
        module procedure globalize_sum_dp, globalize_sum_int
        module procedure globalize_sum_v_dp
      end interface
      interface broadcast
        module procedure broadcast_dp, broadcast_int, broadcast_logical
        module procedure broadcast_sp, broadcast_char
        module procedure broadcast_v_dp
      end interface

      CONTAINS

!------------------------------------------------------------
      subroutine broadcast_dp(scalar)
      real(dp), intent(inout) :: scalar

#ifdef MPI
#ifdef NODAT
           call MPI_Bcast(scalar,1,MPI_double_precision,0,
     $            MPI_Comm_World,MPIerror)
#else
           call MPI_Bcast(scalar,1,DAT_double,0,
     $            MPI_Comm_World,MPIerror)
#endif
#endif
       end subroutine broadcast_dp

!------------------------------------------------------------
      subroutine broadcast_v_dp(a)
      real(dp), dimension(:), intent(inout) :: a

#ifdef MPI
#ifdef NODAT
           call MPI_Bcast(a,size(a),MPI_double_precision,0,
     $            MPI_Comm_World,MPIerror)
#else
           call MPI_Bcast(a,size(a),DAT_double,0,
     $            MPI_Comm_World,MPIerror)
#endif
#endif
       end subroutine broadcast_v_dp

!---------------------------------------------------
      subroutine broadcast_sp(scalar)
      real(sp), intent(inout) :: scalar

#ifdef MPI
#ifdef NODAT
           call MPI_Bcast(scalar,1,MPI_real,0,
     $            MPI_Comm_World,MPIerror)
#else
           call MPI_Bcast(scalar,1,DAT_single,0,
     $            MPI_Comm_World,MPIerror)
#endif
#endif
       end subroutine broadcast_sp

!---------------------------------------------------
      subroutine broadcast_int(scalar)
      integer, intent(inout) :: scalar

#ifdef MPI
           call MPI_Bcast(scalar,1,MPI_Integer,0,
     $            MPI_Comm_World,MPIerror)
#endif
       end subroutine broadcast_int
!---------------------------------------------------
      subroutine broadcast_char(str)
      character(len=*), intent(inout) :: str

#ifdef MPI
           call MPI_Bcast(str,len(str),MPI_Character,0,
     $            MPI_Comm_World,MPIerror)
#endif
       end subroutine broadcast_char
!----------------------------------------------------
      subroutine broadcast_logical(scalar)
      logical, intent(inout) :: scalar

#ifdef MPI
           call MPI_Bcast(scalar,1,MPI_Logical,0,
     $            MPI_Comm_World,MPIerror)
#endif
       end subroutine broadcast_logical

!---------------------------------------------------
!---------------------------------------------------
      subroutine Globalize_sum_dp(local,global)
      real(dp), intent(in) :: local
      real(dp), intent(out) :: global

#ifndef MPI
           global = local
#else
#ifdef NODAT
           call MPI_AllReduce(local,global,1,MPI_double_precision,
     $            MPI_sum,MPI_Comm_World,MPIerror)
#else
           call MPI_AllReduce(local,global,1,DAT_double,
     $            MPI_sum,MPI_Comm_World,MPIerror)
#endif
#endif
       end subroutine Globalize_sum_dp

      subroutine Globalize_sum_v_dp(local,global)
      real(dp), intent(in), dimension(:)  :: local
      real(dp), intent(out), dimension(:) :: global

      integer :: n

      n = size(local)
      if ( n /= size(global))
     $     call die("Globalize_sum_v_dp error")
#ifndef MPI
           global = local
#else
#ifdef NODAT
           call MPI_AllReduce(local,global,n,MPI_double_precision,
     $            MPI_sum,MPI_Comm_World,MPIerror)
#else
           call MPI_AllReduce(local,global,n,DAT_double,
     $            MPI_sum,MPI_Comm_World,MPIerror)
#endif
#endif
       end subroutine Globalize_sum_v_dp

      subroutine Globalize_sum_int(local,global)
      integer, intent(in) :: local
      integer, intent(out) :: global

#ifndef MPI
           global = local
#else
           call MPI_AllReduce(local,global,1,MPI_integer,
     $            MPI_sum,MPI_Comm_World,MPIerror)
#endif
       end subroutine Globalize_sum_int
!-------------------------------------------------------
      subroutine Globalize_max_dp(local,global)
      real(dp), intent(in) :: local
      real(dp), intent(out) :: global

#ifndef MPI
           global = local
#else
#ifdef NODAT
           call MPI_AllReduce(local,global,1,MPI_double_precision,
     $            MPI_max,MPI_Comm_World,MPIerror)
#else
           call MPI_AllReduce(local,global,1,DAT_double,
     $            MPI_max,MPI_Comm_World,MPIerror)
#endif
#endif
       end subroutine Globalize_max_dp

      subroutine Globalize_max_int(local,global)
      integer, intent(in) :: local
      integer, intent(out) :: global

#ifndef MPI
           global = local
#else
           call MPI_AllReduce(local,global,1,MPI_integer,
     $            MPI_max,MPI_Comm_World,MPIerror)
#endif
       end subroutine Globalize_max_int

      end module m_mpi_utils
