! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine diagg( nspin, nuo, maxuo, maxnh, maxnd, 
     &                  maxo, numh, listhptr, listh, numd, 
     &                  listdptr, listd, H, S,
     &                  getD, getPSI, fixspin, qtot, qs, temp, e1, e2,
     &                  eo, qo, Dnew, Enew, ef, efs, Entropy,
     &                  Haux, Saux, psi, nuotot, occtol,
     &                  iscf, neigwanted )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). Gamma-point version.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals local to node
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               ofdensity matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C logical getPSI              : Find and print wavefunctions?
C logical fixspin             : Fix the spin of the system?
C real*8  qtot                : Number of electrons in unit cell
C real*8  qs(nspin)           : Number of electrons in unit cell for each
C                               spin component (if fixed spin option is used)
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C integer iscf                : SCF cycle number
C real*8  occtol              : Occupancy threshold for DM build
C integer neigwanted          : Number of eigenvalues wanted
C *************************** OUTPUT **********************************
C real*8 eo(maxo,nspn)        : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,nspn)        : Occupations of eigenstates
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix
C real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 efs(nspin)           : Fermi energy for each spin
C                               (for fixed spin calculations)
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(nuotot,nuo)     : Auxiliary space for the hamiltonian matrix
C real*8 Saux(nuotot,nuo)     : Auxiliary space for the overlap matrix
C real*8 psi(nuotot,maxuo,nspin) : Auxiliary space for the eigenvectors
C real*8 aux(nuotot)          : Extra auxiliary space
C *************************** UNITS ***********************************
C eo, Enew and ef returned in the units of H.
C *************************** PARALLEL ********************************
C The auxiliary arrays are now no longer symmetry and so the order
C of referencing has been changed in several places to reflect this.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use writewave,     only : writew
      use m_fermid,      only : fermid, fermispin, stepf
      use alloc

#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer 
     &  MPIerror
#endif

      integer
     &  iscf, maxnd, maxnh, maxuo, maxo, nuo, nspin, nuotot,
     &  neigwanted

      integer 
     &  listh(maxnh), numh(nuo), listhptr(nuo),
     &  listd(maxnd), numd(nuo), listdptr(nuo)

      real(dp)
     &  Dnew(maxnd,nspin), e1, e2, ef, Enew(maxnd,nspin), 
     &  Entropy, eo(maxo,nspin), H(maxnh,nspin), qo(maxo,nspin), 
     &  qtot, qs(nspin), S(maxnh), temp, efs(nspin), occtol
     
      real(dp)
     &  Haux(nuotot,nuo), Saux(nuotot,nuo), psi(nuotot,maxuo,nspin) 

      real(dp), dimension(1), parameter :: wk = (/ 1.0_dp /)
      integer, parameter                :: nk = 1

      logical
     &  getD, getPSI, fixspin

      external rdiag

C  Internal variables .............................................
      integer           ie, io, iio, ispin, ix, j, jo, BNode, iie, ind,
     &                  ierror, nd
      real(dp)          ee, eei, qe, qei, rt, t, k(3)

      integer           :: maxnuo, mm
      integer, pointer  :: nuo_LOC(:)
      real(dp), pointer :: psi_tmp(:), paux(:)

C Solve eigenvalue problem
      do ispin = 1,nspin
        call timer( 'r-eigvec', 1 )
        call timer( 'r-buildHS', 1 )
        do io = 1,nuo
          do jo = 1,nuotot
            Saux(jo,io) = 0.0d0
            Haux(jo,io) = 0.0d0
          enddo
        enddo
        do io = 1,nuo
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            Saux(jo,io) = Saux(jo,io) + S(ind)
            Haux(jo,io) = Haux(jo,io) + H(ind,ispin)
          enddo
        enddo
        call timer( 'r-buildHS', 2 )
        call rdiag( Haux, Saux, nuotot, nuo, nuotot, eo(1,ispin), 
     &              psi(1,1,ispin), neigwanted, iscf, ierror )
        call timer( 'r-eigvec', 2 )


C Check error flag and take appropriate action
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
          do io = 1,nuo
            do jo = 1,nuotot
              Saux(jo,io) = 0.0d0
              Haux(jo,io) = 0.0d0
            enddo
          enddo
          do io = 1,nuo
            do j = 1,numh(io)
              ind = listhptr(io) + j
              jo = listh(ind)
              Saux(jo,io) = Saux(jo,io) + S(ind)
              Haux(jo,io) = Haux(jo,io) + H(ind,ispin)
            enddo
          enddo
          call rdiag( Haux, Saux, nuotot, nuo, nuotot, eo(1,ispin), 
     &                psi(1,1,ispin), nuotot, iscf, ierror )
        endif

        if (getPSI) then
          do ix = 1,3
            k(ix) = 0.0d0
          enddo
          call writew(nuotot,nuo,1,k,ispin,
     &                eo(1,ispin),psi(1,1,ispin),.true.)
        endif

      enddo
C Check if we are done ................................................
      if (.not.getD) return

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, maxo, neigwanted, 
     &                  eo, temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, nk, wk, maxo, neigwanted, eo, 
     &               temp, qtot, qo, ef, Entropy )
      endif

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        rt = 1.0d0/t
        do ispin = 1,nspin
          do io = 1,nuotot
            qo(io,ispin) = ( stepf((eo(io,ispin)-e2)*rt) -
     &                       stepf((eo(io,ispin)-e1)*rt)) * 
     &                       2.0d0/dble(nspin)
          enddo
        enddo
      endif
      
      if (nuo.gt.0) then
C New density and energy-density matrices of unit-cell orbitals .......
        nd = listdptr(nuo) + numd(nuo)
        Dnew(1:nd,1:nspin) = 0.0d0
        Enew(1:nd,1:nspin) = 0.0d0
      endif

      call timer( 'r-buildD', 1 )
C Global operation to form new density matrix
      nullify( nuo_LOC )
      call re_alloc( nuo_LOC, 0, Nodes-1, 'nuo_LOC', 'diagg' )
#ifdef MPI
      call MPI_Allgather( nuo, 1, MPI_INTEGER, nuo_LOC, 1,
     &                    MPI_INTEGER, MPI_Comm_World, MPIerror )
#else
      nuo_LOC(0) = nuo
#endif
      maxnuo = 0
      do BNode=0, Nodes-1
        maxnuo = max(nuo_LOC(BNode),maxnuo)
      enddo
      nullify( PSI_TMP )
      call re_alloc( PSI_TMP, 1, nuotot*maxnuo*nspin, 'PSI_TMP',
     &               'diagg' )

      do Bnode=0, Nodes-1
        if (BNode.eq.Node) then
          mm = 0
          do ispin= 1, nspin
            do io= 1, nuo
              PSI_TMP(mm+1:mm+nuotot) = PSI(1:nuotot,io, ispin)
              mm = mm + nuotot
            enddo
          enddo
        endif
#ifdef MPI
        call MPI_Bcast( PSI_TMP, nuotot*nuo_LOC(Bnode)*nspin,
     &                  MPI_double_precision, BNode, MPI_Comm_World,
     &                  MPIerror )
#endif
        do ispin = 1,nspin
          do io = 1,nuo
            call LocalToGlobalOrb(io,Node,Nodes,iio)
            mm   = (ispin-1)*nuo_LOC(Bnode)*nuotot
            do iie = 1,nuo_LOC(Bnode)
              call  LocalToGlobalOrb(iie,BNode,Nodes,ie)
              qe = qo(ie,ispin)
              if (abs(qe).gt.occtol) then
                paux => PSI_TMP(mm+1:mm+nuotot)
                ee  = qe*eo(ie,ispin)
                qei = qe*paux(iio)
                eei = ee*paux(iio)
                do j = 1,numd(io)
                  ind = listdptr(io) + j
                  jo  = listd(ind)
                  Dnew(ind,ispin) = Dnew(ind,ispin) + qei*paux(jo)
                  Enew(ind,ispin) = Enew(ind,ispin) + eei*paux(jo)
                enddo
              endif
              mm = mm + nuotot
            enddo
          enddo
        enddo
      enddo

      call timer( 'r-buildD', 2 )
      call de_alloc( nuo_LOC, 'nuo_LOC', 'diagg' )
      call de_alloc( PSI_TMP, 'PSI_TMP', 'diagg' )
      end
