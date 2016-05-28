MODULE MPI_INSTR
!
! This is an interface to supplant some MPI routines
! in order to time-profile them. 
! J.M.Soler. May.2009
! A. Garcia, June 2015
!
  USE MPI_INTERFACES

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
        public :: mpi_finalize


  private :: timer_mpi

CONTAINS

        subroutine timer_mpi(prog,iopt)
          character(len=*), intent(in) :: prog
          integer, intent(in)          :: iopt
          
          include 'timer_mpi_handler.inc'
        end subroutine timer_mpi

END MODULE MPI_INSTR
