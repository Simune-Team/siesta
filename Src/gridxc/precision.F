!!@LICENSE
!
      module precision
!
!     Precision handling
!
#ifdef MPI
      use mpi, only: MPI_real, MPI_double_precision
#endif

      implicit none
      
      integer, parameter :: i8b = selected_int_kind(18)

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
!
!     Select precision of certain parts of the program
!     These are set through preprocessor directives.
!
#ifdef MPI
#ifdef GRID_SP
      integer, parameter :: grid_p = sp
      integer, parameter :: MPI_grid_real = MPI_real
#elif defined(GRID_DP)
      integer, parameter :: grid_p = dp
      integer, parameter :: MPI_grid_real = MPI_double_precision
#else
      integer, parameter :: grid_p = dp
      integer, parameter :: MPI_grid_real = MPI_double_precision
#endif

#else   

#ifdef GRID_SP
      integer, parameter :: grid_p = sp
#elif defined(GRID_DP)
      integer, parameter :: grid_p = dp
#else
      integer, parameter :: grid_p = dp
#endif

#endif
      end module precision

