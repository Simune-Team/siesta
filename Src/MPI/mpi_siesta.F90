MODULE MPI_SIESTA

! For this PEXSI version, temporarily removed timing versions of some
! MPI routines.

#undef MPI

! The following construction allows to supplant MPI_Comm_World within SIESTA,
! and to use it as a subroutine with its own internal MPI communicator.
! JMS. Oct.2010, AG, March 2013

  USE MPI, true_MPI_Comm_World => MPI_Comm_World
  integer, public :: MPI_Comm_World = true_MPI_Comm_World
  public :: true_MPI_Comm_World


#ifdef GRID_SP
        integer, parameter :: MPI_grid_real   = MPI_real
#elif defined(GRID_DP)
        integer, parameter :: MPI_grid_real   = MPI_double_precision
#else
        integer, parameter :: MPI_grid_real   = MPI_double_precision
#endif


!
!   Export explicitly some symbols to help some versions of
!   the PGI compiler. It does not consider them public by default
!
        public :: mpi_real
        public :: mpi_complex
        public :: mpi_double_complex
        public :: mpi_double_precision
        public :: mpi_2double_precision
        public :: mpi_integer, mpi_character, mpi_logical
        public :: mpi_integer8
        public :: mpi_maxloc, mpi_sum, mpi_max, mpi_lor
        public :: mpi_status_size
        public :: mpi_grid_real
        public :: mpi_group_null, mpi_comm_null, mpi_proc_null

END MODULE MPI_SIESTA
