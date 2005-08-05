      subroutine savepsi(psiprev, psi, maxpsi, nuo, nuotot, nocc)
C *********************************************************************
C Copies psi into psiprev allowing for the fact that the parallel form
C requires a re-distribution of the data.
C Written by J.D. Gale, March 2000
C **************************** INPUT ********************************** 
C real*8  psi(maxpsi)         : Wavefunctions in current k point
C real*8  maxpsi              : Array dimensions for psi and psiprev
C integer nuo                 : Number of basis orbitals in the unit cell
C integer nuotot              : Third dimension of xij
C integer nocc                : number of occupied states
C *************************** OUTPUT **********************************
C real*8  psiprev(maxpsi)     : Wavefunctions in previous k point
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GetNodeOrbs, LocalToGlobalOrb
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
 
      integer nuo, nuotot, maxpsi, nocc

      real(dp)
     .  psiprev(2,nuo,nuotot), 
     .  psi(2,nuotot,nuo)

C**** Internal variables ***********************************************

      integer 
     .  iuo, juo
#ifdef MPI
      integer 
     .  MPIerror, iuog, juog, n, noccloc
      real(dp), dimension(:,:,:), allocatable, save ::
     .  psitmp
#endif

#ifdef MPI
      do n = 0,Nodes-1

C Broadcast copy of psi on node n to all other nodes
        call GetNodeOrbs(nocc,n,Nodes,noccloc)
        allocate(psitmp(2,nuotot,noccloc))
        call memory('A','D',2*nuotot*noccloc,'detover')
        if (Node .eq. n) then
          do iuo = 1,noccloc
            do juo = 1,nuotot
              psitmp(1,juo,iuo) = psi(1,juo,iuo)
              psitmp(2,juo,iuo) = psi(2,juo,iuo)
            enddo
          enddo
        endif
        call MPI_Bcast(psitmp(1,1,1),2*nuotot*noccloc,
     .    MPI_double_precision,n,MPI_Comm_World,MPIerror)

C Save local part of psiprev
        do iuo = 1,noccloc
          call LocalToGlobalOrb(iuo,n,Nodes,iuog)
          do juo = 1,nuo
            call LocalToGlobalOrb(juo,Node,Nodes,juog)
            psiprev(1,juo,iuog) = psitmp(1,juog,iuo)
            psiprev(2,juo,iuog) = psitmp(2,juog,iuo)
          enddo
        enddo

        call memory('D','D',size(psitmp),'detover')
        deallocate(psitmp)

      enddo
#else
C Straight serial copy
      do iuo = 1,nocc
        do juo = 1,nuotot
          psiprev(1,juo,iuo) = psi(1,juo,iuo)
          psiprev(2,juo,iuo) = psi(2,juo,iuo)
        enddo
      enddo
#endif

      end
