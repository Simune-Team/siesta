MODULE MPI_SIESTA
!
! This is an interface to supplant some MPI routines called by siesta,
! in order to time-profile them. J.M.Soler. May.2009
!
  USE MPI_INTERFACES,  &! Previously called MPI_SIESTA
    trueMPI_BARRIER    => MPI_BARRIER,    & ! Renamed to avoid conflicts
    trueMPI_COMM_RANK  => MPI_COMM_RANK,  &
    trueMPI_COMM_SIZE  => MPI_COMM_SIZE,  &
    trueMPI_COMM_SPLIT => MPI_COMM_SPLIT, &
    trueMPI_GET_COUNT  => MPI_GET_COUNT,  &
    trueMPI_INIT       => MPI_INIT,       &
    trueMPI_WAIT       => MPI_WAIT,       &
    trueMPI_WAITALL    => MPI_WAITALL

  USE TIMER_MPI_M, only: timer_mpi

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
        public :: mpi_comm_world
        public :: mpi_comm_self
        public :: mpi_grid_real
        public :: mpi_finalize


  PUBLIC :: MPI_BARRIER
  INTERFACE MPI_BARRIER
    MODULE PROCEDURE myMPI_BARRIER
  END INTERFACE

  PUBLIC :: MPI_COMM_RANK
  INTERFACE MPI_COMM_RANK
    MODULE PROCEDURE myMPI_COMM_RANK
  END INTERFACE

  PUBLIC :: MPI_COMM_SIZE
  INTERFACE MPI_COMM_SIZE
    MODULE PROCEDURE myMPI_COMM_SIZE
  END INTERFACE

  PUBLIC :: MPI_COMM_SPLIT
  INTERFACE MPI_COMM_SPLIT
    MODULE PROCEDURE myMPI_COMM_SPLIT
  END INTERFACE

  PUBLIC :: MPI_GET_COUNT
  INTERFACE MPI_GET_COUNT
    MODULE PROCEDURE myMPI_GET_COUNT
  END INTERFACE

  PUBLIC :: MPI_INIT
  INTERFACE MPI_INIT
    MODULE PROCEDURE myMPI_INIT
  END INTERFACE

  PUBLIC :: MPI_WAIT
  INTERFACE MPI_WAIT
    MODULE PROCEDURE myMPI_WAIT
  END INTERFACE

  PUBLIC :: MPI_WAITALL
  INTERFACE MPI_WAITALL
    MODULE PROCEDURE myMPI_WAITALL
  END INTERFACE

CONTAINS

  SUBROUTINE myMPI_BARRIER(COMM, IERROR) 
    INTEGER, INTENT(IN)  :: COMM
    INTEGER, INTENT(OUT) :: IERROR 
    external MPI_BARRIER
    call timer_mpi('MPI_BARRIER',1)
    call MPI_BARRIER(COMM, IERROR)
    call timer_mpi('MPI_BARRIER',2)
  END SUBROUTINE myMPI_BARRIER
          
  SUBROUTINE myMPI_COMM_RANK(COMM, RANK, IERROR)
    INTEGER, INTENT(IN)  :: COMM
    INTEGER, INTENT(OUT) :: RANK
    INTEGER, INTENT(OUT) :: IERROR 
    external MPI_COMM_RANK
    call timer_mpi('MPI_COMM_RANK',1)
    call MPI_COMM_RANK(COMM, RANK, IERROR)
    call timer_mpi('MPI_COMM_RANK',2)
  END SUBROUTINE myMPI_COMM_RANK
          
  SUBROUTINE myMPI_COMM_SIZE(COMM, SIZE, IERROR)
    INTEGER, INTENT(IN)  :: COMM
    INTEGER, INTENT(OUT) :: SIZE
    INTEGER, INTENT(OUT) :: IERROR 
    external MPI_COMM_SIZE
    call timer_mpi('MPI_COMM_SIZE',1)
    call MPI_COMM_SIZE(COMM, SIZE, IERROR)
    call timer_mpi('MPI_COMM_SIZE',2)
  END SUBROUTINE myMPI_COMM_SIZE
          
  SUBROUTINE myMPI_COMM_SPLIT(COMM, COLOR, KEY, NEWCOMM, IERROR)
    INTEGER, INTENT(IN)  :: COMM
    INTEGER, INTENT(IN)  :: COLOR
    INTEGER, INTENT(IN)  :: KEY
    INTEGER, INTENT(OUT) :: NEWCOMM
    INTEGER, INTENT(OUT) :: IERROR 
    external MPI_COMM_SPLIT
    call timer_mpi('MPI_COMM_SPLIT',1)
    call MPI_COMM_SPLIT(COMM, COLOR, KEY, NEWCOMM, IERROR)
    call timer_mpi('MPI_COMM_SPLIT',2)
  END SUBROUTINE myMPI_COMM_SPLIT
          
  SUBROUTINE myMPI_GET_COUNT(STATUS, DATATYPE, COUNT, IERROR)
    USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
    INTEGER, INTENT(IN)  :: STATUS(MPI_STATUS_SIZE)
    INTEGER, INTENT(IN)  :: DATATYPE
    INTEGER, INTENT(OUT) :: COUNT
    INTEGER, INTENT(OUT) :: IERROR
    external MPI_GET_COUNT
    call timer_mpi('MPI_GET_COUNT',1)
    call MPI_GET_COUNT(STATUS, DATATYPE, COUNT, IERROR)
    call timer_mpi('MPI_GET_COUNT',2)
  END SUBROUTINE myMPI_GET_COUNT
          
  SUBROUTINE myMPI_INIT(IERROR)
    INTEGER, INTENT(OUT) :: IERROR 
    external MPI_INIT
    call timer_mpi('MPI_INIT',1)
    call MPI_INIT(IERROR)
    call timer_mpi('MPI_INIT',2)
  END SUBROUTINE myMPI_INIT
          
  SUBROUTINE myMPI_WAIT(REQUEST, STATUS, IERROR)
    USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
    INTEGER, INTENT(INOUT) :: REQUEST
    INTEGER, INTENT(OUT) :: STATUS(MPI_STATUS_SIZE)
    INTEGER, INTENT(OUT) :: IERROR 
    external MPI_WAIT
    call timer_mpi('MPI_WAIT',1)
    call MPI_WAIT(REQUEST, STATUS, IERROR)
    call timer_mpi('MPI_WAIT',2)
  END SUBROUTINE myMPI_WAIT
          
  SUBROUTINE myMPI_WAITALL( &
    COUNT, ARRAY_OF_REQUESTS, ARRAY_OF_STATUSES, IERROR)
    USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
    INTEGER, INTENT(IN)  :: COUNT
    INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
    INTEGER, INTENT(OUT) :: ARRAY_OF_STATUSES(MPI_STATUS_SIZE,*)
    INTEGER, INTENT(OUT) :: IERROR 
    external MPI_WAITALL
    call timer_mpi('MPI_WAITALL',1)
    call MPI_WAITALL(COUNT, ARRAY_OF_REQUESTS, ARRAY_OF_STATUSES, IERROR)
    call timer_mpi('MPI_WAITALL',2)
  END SUBROUTINE myMPI_WAITALL
          
END MODULE MPI_SIESTA
