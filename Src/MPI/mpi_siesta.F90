MODULE MPI_SIESTA

#undef MPI

  USE MPI, true_MPI_Comm_World => MPI_Comm_World

#ifdef GRID_SP
        integer, parameter :: MPI_grid_real   = MPI_real
#elif defined(GRID_DP)
        integer, parameter :: MPI_grid_real   = MPI_double_precision
#else
        integer, parameter :: MPI_grid_real   = MPI_double_precision
#endif

  integer, public :: MPI_Comm_World = true_MPI_Comm_World
  public :: true_MPI_Comm_World
!
!   Export explicitly some symbols to help some versions of
!   the PGI compiler, which do not consider them public by default
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

END MODULE MPI_SIESTA
