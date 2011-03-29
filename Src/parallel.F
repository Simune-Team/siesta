!!@LICENSE
!
      module parallel

C  Parallelisation related global parameters
C
C  integer BlockSize = this is the blocking factor used to divide up
C                      the arrays over the processes for the Scalapack
C                      routines. Setting this value is a compromise
C                      between dividing up the orbitals over the processors
C                      evenly to achieve load balancing and making the
C                      local work efficient. Typically a value of about
C                      10 is good, but optimisation may be worthwhile.
C                      A value of 1 is very bad for any number of processors
C                      and a large value may also be less than ideal.
C
C  integer ProcessorY = second dimension of processor grid in mesh point
C                       parallelisation - note that the first dimension
C                       is determined by the total number of processors
C                       in the current job. Also note that this number
C                       must be a factor of the total number of processors.
C                       Furthermore on many parallel machines (e.g. T3E)
C                       this number must also be a power of 2.
C
      implicit none

      integer, save :: Node = 0
      integer, save :: Nodes = 1
      integer, save :: BlockSize  = 24
      integer, save :: ProcessorY = 1
      
      logical, save :: IOnode = .true.
      logical, save :: ParallelOverK = .false.
      logical, save :: ResetFirstCall = .false.

      public

      contains

      subroutine parallel_init()  ! Initializes Node, Nodes, and IOnode
#ifdef MPI
      use mpi_siesta, only: MPI_Comm_World
      logical,save:: initialized=.false.
      integer :: MPIerror
      if (.not.initialized) then
        call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
        call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
        IOnode = (Node==0)
        initialized = .true.
      end if
#endif
      end subroutine parallel_init

      end module parallel
!
!     convenience function
! 
      function io_node()
      use parallel, only: node
      logical :: io_node
      
      io_node = ( node == 0)
      end function io_node
