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
      subroutine reordpsi(psiprev, psi, nuo, nuotot, nocc)
C *********************************************************************
C Reorders the coefficient of the bands, so the bands whose bands index
C range between 1 and numincbands are the bands included for
C wannierization.
C Written by J. Junquera, October 2013, based on subroutine savepsi,
C by J.D. Gale, March 2000
C **************************** INPUT ********************************** 
C real*8  psi(2,nuotot,nuo)     : Wavefunctions in current k point
C integer nuo                   : Number of (local) orbitals in the unit cell
C integer nuotot                : Number of orbitals in the unit cell
C integer nocc                  : number of occupied states
C *************************** OUTPUT **********************************
C real*8  psiprev(2,nuo,nuotot) : Wavefunctions in previous k point
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
      use precision
      use parallel,           only : Node, Nodes
      use parallelsubs,       only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,       only : GlobalToLocalOrb
      use alloc,              only : re_alloc, de_alloc
      use m_siesta2wannier90, only : isexcluded
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
 
      integer nuo, nuotot, nocc, jband, jbandloc

      real(dp)
     .  psiprev(2,nuotot,nuo), 
     .  psi(2,nuotot,nuo)

C**** Internal variables ***********************************************

      integer 
     .  iuo, juo, iband
#ifdef MPI
      integer 
     .  MPIerror, iuog, juog, n, noccloc
      real(dp), dimension(:,:,:), pointer ::  psitmp
#endif

#ifdef MPI
! AG
! Allocate as in detover, using the number of orbitals on the first node,
! as some of the nodes might have zero orbitals!
      call GetNodeOrbs(nocc,0,Nodes,noccloc)
      nullify( psitmp )
      call re_alloc( psitmp, 1, 2, 1, nuotot, 1, noccloc,
     &                 name='psitmp', routine='reordpsi' )

      jband = 0
      do n = 0,Nodes-1

C Broadcast copy of psi on node n to all other nodes
! Note how not all the bands are broadcasted. 
! Only the "occupied" bands in the corresponding node. 
        call GetNodeOrbs(nocc,n,Nodes,noccloc)

        if (Node .eq. n) then
           psitmp(1:2,1:nuotot,1:noccloc) = psi(1:2,1:nuotot,1:noccloc)
        endif
        call MPI_Bcast(psitmp(1,1,1),2*nuotot*noccloc,
     .    MPI_double_precision,n,MPI_Comm_World,MPIerror)

C Save local part of psiprev
        do iuo = 1,noccloc
          call LocalToGlobalOrb(iuo,n,Nodes,iuog)
!         Select here if the local occupied band is included or not.
!         If it is included, then a new index (jband) that ranges
!         between 1 and numincbands(ispin) is assigned 
          if( .not.  isexcluded(iuog) ) then
            jband = jband + 1
            call GlobalToLocalOrb(jband,Node,Nodes,jbandloc)
            if( jbandloc .ne. 0) then
!! For debugging
!              write(6,*)'Node = ', Node, n, iuo, iuog, jband, jbandloc
!! End debugging
              do juo = 1,nuotot
                psiprev(1,juo,jbandloc) = psitmp(1,juo,iuo)
                psiprev(2,juo,jbandloc) = psitmp(2,juo,iuo)
              enddo
            endif
          endif
        enddo


      enddo
      call de_alloc( psitmp,  name='psitmp' )

#else
C Straight serial copy
      jband = 0
      do iband = 1, nocc
        if (.not. isexcluded(iband)) then
          jband = jband + 1
          psiprev(1,:,jband) = psi(1,:,iband)
          psiprev(2,:,jband) = psi(2,:,iband)
        endif
      enddo
#endif

      end
