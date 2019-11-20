! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! Stand-alone 'bye' routine for use by libraries and
! low-level modules.
!
! Each program using the module or library needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.

!---------------------------------------------------------
      subroutine bye(str)

#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

      integer  :: Node
      external :: pxfflush
#ifdef MPI
      integer rc, MPIerror
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
         call pxfflush(6)
      endif

#ifdef MPI
      call MPI_Finalize(rc)
#endif
      stop
      
      end subroutine bye
