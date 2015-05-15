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
      subroutine reordpsi(psisave, psi, nuo, nuotot, nocc, ninc)
C *********************************************************************
C Reorders the coefficient of the bands, so the bands whose bands index
C range between 1 and numincbands are the bands included for
C wannierization.
C Written by J. Junquera, October 2013, based on subroutine savepsi,
C by J.D. Gale, March 2000
C **************************** INPUT ********************************** 
C real*8  psi(2,nuotot,nuo)      : Wavefunctions in current k point
C integer nuo                    : Number of (local) orbitals in the unit cell
C integer nuotot                 : Number of orbitals in the unit cell
C integer nocc                   : Number of bands considered for
C                                  Wannierization (before excluding bands)
C integer ninc                   : Number of maximum bands whose
C                                  coefficients will be stored in a node
C *************************** OUTPUT **********************************
C real*8  psisave(2,nuotot,ninc) : Wavefunctions saved in the matrix of
C                                  coefficients
C                                  Both in psi and psisave,
C                                  First  index: real or complex value
C                                  Second index: atomic orbital
C                                  Third  index: band index
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
      use precision
      use m_siesta2wannier90, only : isexcluded
#ifdef MPI
      use alloc,              only : re_alloc, de_alloc
      use parallel,           only : Node, Nodes, BlockSize
      use parallelsubs,       only : GetNodeOrbs, LocalToGlobalOrb
      use parallelsubs,       only : GlobalToLocalOrb
      use mpi_siesta
      use m_orderbands,       only : index_included_band
      use m_orderbands,       only : node_included_band
      use m_orderbands,       only : band_index_in_node
      use m_orderbands,       only : which_band_in_node
#endif
!     For debugging
      use sys,                only : die
!     End debugging
      implicit none
 
      integer nuo, nuotot, nocc, ninc, jband, jbandloc

      real(dp)
     .  psisave(2,nuotot,ninc), 
     .  psi(2,nuotot,nuo)

C**** Internal variables ***********************************************

      integer 
     .  iuo, juo, iband
#ifdef MPI
      integer 
     .  MPIerror, iuog, juog, n, noccloc
      real(dp), dimension(:,:,:), pointer ::  psitmp
#endif

!! For debugging
!#ifdef MPI
!      write(6,'(a,3i5)')'reordpsi: Node, nocc, ninc', Node, nocc, ninc
!#endif
!! End debugging
!
#ifdef MPI
! AG
! Allocate as in detover, using the number of orbitals on the first node,
! as some of the nodes might have zero orbitals!
      call GetNodeOrbs(nocc,0,Nodes,noccloc)

      nullify( psitmp )
      call re_alloc( psitmp, 1, 2, 1, nuotot, 1, noccloc,
     &                 name='psitmp', routine='reordpsi' )

!! For debugging
!      if( Node .eq. 4 ) then
!        do iuo = 1, ninc
!          write(6,'(a,4i5)')
!     .      ' reordpsi: Node, iuo, index_included, node_included =', 
!     .      Node, iuo, index_included_band(iuo), node_included_band(iuo)
!        enddo 
!      endif
!
!      do jband = 1, ninc
!        write(6,'(a,3i5)')' Node, index_in_node, orbital = ',
!     .    Node, jband, which_band_in_node(Node,jband)
!      enddo
!! End debugging

      jband = 0
      do n = 0, Nodes-1

!       Broadcast copy of psi on node n to all other nodes
!       Note how not all the bands are broadcasted. 
!       Only the "occupied" bands in the corresponding node. 
        call GetNodeOrbs(nocc,n,Nodes,noccloc)

!!       For debugging
!           write(6,'(a,4i5)')' reordpsi: Node, n, nocc, noccloc = ', 
!     .                                   Node, n, nocc, noccloc
!!       End debugging

        if (Node .eq. n) then
           psitmp(1:2,1:nuotot,1:noccloc) = psi(1:2,1:nuotot,1:noccloc)
        endif
        call MPI_Bcast(psitmp(1,1,1),2*nuotot*noccloc,
     .    MPI_double_precision,n,MPI_Comm_World,MPIerror)

!       Save local part of psisave
        do iuo = 1, noccloc
          call LocalToGlobalOrb(iuo,n,Nodes,iuog)

!!         For debugging
!            write(6,'(a,4i5)')' reordpsi: Node, n, iuo, iuog = ', 
!     .                                    Node, n, iuo, iuog
!!         End debugging

!         Select here if the local occupied band is included or not.
!         If it is included, then a new index (jband) that ranges
!         between 1 and numincbands(ispin) is assigned 
          if( .not.  isexcluded(iuog) ) then
            jband    = jband + 1
!           The coefficients of the included band will be stored
!           in Node (node_included_band(jband)),
!           and will be the band_index_in_node(jband)-th band
!           stored on that Node
            jbandloc = band_index_in_node(jband)
            if( node_included_band(jband) .eq. Node ) then
!! For debugging
!              write(6,'(a,8i5)')
!     .          ' reordpsi: Node, n, iuo, iuog, jband = ', 
!     .          Node, n, iuo, iuog, jband, jbandloc,
!     .          index_included_band(jband), node_included_band(jband)
!! End debugging
              do juo = 1, nuotot
                psisave(1,juo,jbandloc) = psitmp(1,juo,iuo)
                psisave(2,juo,jbandloc) = psitmp(2,juo,iuo)
              enddo ! atomic orbitals
            endif   ! If the local Node is the same where the included
                    !    band should be stored
          endif     ! Is the band included for wannierization?
        enddo       ! local bands on a given node
      enddo         ! nodes
      call de_alloc( psitmp,  name='psitmp' )
#else
!     Straight serial copy
      jband = 0
      do iband = 1, nocc
        if (.not. isexcluded(iband)) then
          jband = jband + 1
          psisave(1,:,jband) = psi(1,:,iband)
          psisave(2,:,jband) = psi(2,:,iband)
        endif
      enddo
#endif

!!     For debugging
!      call die()
!!     End debugging

      end subroutine reordpsi

