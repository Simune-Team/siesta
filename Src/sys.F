      module sys
!
!     Termination and messaging routines, MPI aware
!
      public :: die, bye, message

      CONTAINS

!--------------------------------------------------------
      subroutine message(str)

      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str


      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
      endif

      end subroutine message
!------------------------------------------------------------------
!--------------------------------------------------------
      subroutine die(str)

      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

#ifdef MPI
      integer MPIerror
#endif

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a,i4)') 'Stopping Program from Node: ', Node
!!                                       endif
#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      stop
#endif
      end subroutine die
!------------------------------------------------------------------
      subroutine bye(str)

      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

#ifdef MPI
      integer rc
#endif


! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Requested End of Run. Bye!!'
!!                                       endif
#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine bye
!-----------------------------------------------------------------
      end module sys

