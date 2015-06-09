module gridxc_config
!
!  Contains global data for internal operation
!  of the library.
!  This makes it non-thread-safe for now
!
! MPI variables
#ifdef MPI
  integer,             save :: gridxc_comm
  integer,public,      save :: gridxc_myNode= -1, gridxc_totNodes=-1
  logical,             save :: mpi_initialized = .false.
#else
  integer,public,      save :: gridxc_myNode= 0, gridxc_totNodes=1
#endif
!
! Timing options
!
 logical, public :: gridxc_use_walltime = .false.
 logical, public :: gridxc_time_mpi_calls = .false.

CONTAINS

#ifdef MPI
SUBROUTINE gridxc_mpi_init(comm)
integer, intent(in) :: comm

integer :: mpierr
if (.not. mpi_initialized) then
   gridxc_comm = comm
   call MPI_Comm_Size(comm,gridxc_totNodes,mpierr)
   call MPI_Comm_Rank(comm,gridxc_myNode,mpierr)
   mpi_initialized = .true.
endif
end subroutine gridxc_mpi_init
#endif


SUBROUTINE gridxc_init(comm,            &
                       use_walltime,   &
                       time_mpi_calls)
!                       unit_for_output,...)

integer, intent(in), optional :: comm
!
! internal operation of the library only
logical, intent(in), optional :: use_walltime
logical, intent(in), optional :: time_mpi_calls

#ifdef MPI
if (present(comm)) then
   call gridxc_mpi_init(comm)
else
   call gridxc_mpi_init(MPI_Comm_World)
endif
#else
   gridxc_myNode = 0
   gridxc_totNodes = 1
#endif

gridxc_use_walltime = .false.
if (present(use_walltime)) then
    gridxc_use_walltime = use_walltime
endif

gridxc_time_mpi_calls = .false.
if (present(time_mpi_calls)) then
    gridxc_time_mpi_calls = time_mpi_calls
endif
   
end subroutine gridxc_init

end module gridxc_config

