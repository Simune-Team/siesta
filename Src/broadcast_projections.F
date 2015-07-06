! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine broadcast_projections()
!
!     Globalizes the information on the projection centers for the 
!     maximally localized Wannier functions. 
!
!     Javier Junquera, June 2013, based on 
!     broadcast_basis by Alberto Garcia

      use trialorbitalclass        
      use parallel,           only : Node, Nodes
      use m_siesta2wannier90, only : numproj, projections

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifndef MPI
!
!     Do nothing...
!
      end subroutine broadcast_projections
#else
      integer MPIerror

      integer is, i
      type(trialorbital), pointer  :: proj

#ifdef DEBUG
      call write_debug( '  PRE broadcast_projections' )
#endif

      if ( Nodes .eq. 1 ) return

      call MPI_Bcast(numproj,1,MPI_integer,0,MPI_Comm_World,MPIerror)

      if ( Node .ne. 0 ) allocate (projections(numproj))

      do is = 1, numproj
        proj => projections(is)
        call MPI_Bcast(proj%center,3,MPI_double_precision,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%zaxis,3,MPI_double_precision,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%xaxis,3,MPI_double_precision,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%yaxis,3,MPI_double_precision,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%zovera,1,MPI_double_precision,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%r,1,MPI_integer,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%l,1,MPI_integer,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%mr,1,MPI_integer,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%rcut,1,MPI_double_precision,
     $                 0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(proj%lmax,1,MPI_integer,
     $                 0,MPI_Comm_World,MPIerror)
      enddo

      end subroutine broadcast_projections

#endif



