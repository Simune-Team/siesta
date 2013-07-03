!!@LICENSE

      module sys
!
!     Termination and messaging routines, MPI aware
!
      implicit none

      public :: die      ! Prints an error message and calls MPI_Abort
      public :: bye      ! Prints an error message and calls MPI_Finalize
      public :: message  ! Prints a message string if node==0

      private

      integer, parameter :: Node = 0

      CONTAINS

      subroutine message(str)

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
!
!--------------------------------------------------
      subroutine die(str)

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
            write(0,'(a)') trim(str)
         endif
         write(6,'(a,i4)') 'Stopping Program from Node: ', Node
         write(0,'(a,i4)') 'Stopping Program from Node: ', Node
!!                                       endif
         if (Node .eq. 0) then
            call pxfflush(6)
            call pxfflush(0)
         endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call pxfabort()
#endif
      end subroutine die

!---------------------------------------------------------
      subroutine bye(str)

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

         if (Node .eq. 0) call pxfflush(6)

#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine bye

      end module sys

