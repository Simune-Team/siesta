      module sys

      CONTAINS

      subroutine die(str)
#ifdef MPI
      use mpi
#endif

      character(len=*), intent(in), optional   :: str

      integer Node
#ifdef MPI
      integer MPIerror, rc
#endif

#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Stopping Program'
         write(6,*) 1.0/real(Node)
      endif
#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
!!      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine die

      subroutine bye(str)
#ifdef MPI
      use mpi
#endif

      character(len=*), intent(in), optional   :: str

      integer Node
#ifdef MPI
      integer MPIerror, rc
#endif

#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Requested End of Run. Bye!!'
      endif
#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine bye

      end module sys

