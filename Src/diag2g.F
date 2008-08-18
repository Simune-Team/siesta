! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine diag2g( nuo, no, maxnh, maxnd, maxo,
     .                   numh, listhptr, listh, numd, listdptr,
     .                   listd, H, S, getD, qtot, temp, e1, e2,
     .                   eo, qo, Dnew, Enew, ef, Entropy,
     .                   Haux, Saux, psi, aux, nuotot, occtol,
     .                   iscf, neigwanted )
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for non-colinear spin at gamma point.
C Writen by J.Soler, May and August 1998.
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals on local node
C integer no                  : Number of basis orbitals
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of Dnew and Enew
C integer maxo                : First dimension of eo and qo
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,4)          : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C real*8  qtot                : Number of electrons in unit cell
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C real*8  occtol              : Occupancy threshold for DM build
C integer iscf                : SCF cycle number
C integer neigwanted          : Number of eigenvalues wanted
C *************************** OUTPUT **********************************
C real*8 eo(maxo*2)           : Eigenvalues
C real*8 qo(maxo*2)           : Occupations of eigenstates
C real*8 Dnew(maxnd,4)        : Output Density Matrix
C real*8 Enew(maxnd,4)        : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,2,nuotot,2,nuo): Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,2,nuotot,2,nuo): Auxiliary space for the overlap matrix
C real*8 psi(2,2,nuotot,2*no) : Auxiliary space for the eigenvectors
C real*8 aux(5,2*nuotot)      : Extra auxiliary space
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
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
      use parallel,     only : Node, Nodes, BlockSize
      use parallelsubs, only : LocalToGlobalOrb
      use m_fermid,     only : fermid, stepf
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer
     .  maxo, maxnd, maxnh, no, nuo, nuotot, iscf, neigwanted

      integer 
     .  listh(maxnh), numh(nuo), listhptr(nuo),
     .  listd(maxnd), numd(nuo), listdptr(nuo)

      real(dp)
     .  Dnew(maxnd,4), e1, e2, ef, Enew(maxnd,4), Entropy, eo(maxo*2),
     .  H(maxnh,4), qo(maxo*2), qtot, S(maxnh), temp, occtol
     
      real(dp)
     .  aux(5,2*nuotot), Haux(2,2,nuotot,2,nuo),
     .  psi(2,2,nuotot,2*no), Saux(2,2,nuotot,2,nuo)

      logical
     .  getD

      external cdiag

C  Internal variables .............................................
      integer           BNode, BTest, ie, ierror, iie, iio,
     .                  ind, io, j, jo, nd
      real(dp)          ee, pipj, qe, t
C  ....................

C Initialize Hamiltonian and overlap matrices in full format .....
C Index i is for real/imag parts
C Indices is and js are for spin components
C Indices iuo and juo are for orbital components:
C Haux(i,js,juo,is,iuo) = <js,juo|H|is,iuo>
      Saux = 0.0d0
      Haux = 0.0d0

C Copy S,H matrices from sparse to full format ....................
C Convention: ispin=1 => H11, ispin=2 => H22, ispin=3 => Real(H12),
C             ispin=4 => Imag(H12)
      do io = 1,nuo
        do j = 1,numh(io)
          ind = listhptr(io) + j
          jo = listh(ind)
          Saux(1,1,jo,1,io) =  S(ind)
          Saux(1,2,jo,2,io) =  S(ind)
          Haux(1,1,jo,1,io) =  H(ind,1)
          Haux(1,2,jo,2,io) =  H(ind,2)
          Haux(1,1,jo,2,io) =  H(ind,3)
          Haux(2,1,jo,2,io) = -H(ind,4)
        enddo
      enddo

C Hermiticity
      do io = 1,no
        do jo = 1,no
          Haux(1,2,jo,1,io) =   Haux(1,1,io,2,jo)
          Haux(2,2,jo,1,io) = - Haux(2,1,io,2,jo)
        enddo
      enddo

C Solve the eigenvalue problem .......................................
      call cdiag(Haux,Saux,2*nuotot,2*nuotot,2*nuo,eo,psi,
     .           neigwanted,iscf,ierror)

C Check error flag and take appropriate action
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
        Saux = 0.0d0
        Haux = 0.0d0
        do io = 1,nuo
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            Saux(1,1,jo,1,io) =  S(ind)
            Saux(1,2,jo,2,io) =  S(ind)
            Haux(1,1,jo,1,io) =  H(ind,1)
            Haux(1,2,jo,2,io) =  H(ind,2)
            Haux(1,1,jo,2,io) =  H(ind,3)
            Haux(2,1,jo,2,io) = -H(ind,4)
          enddo
        enddo
        do io = 1,no
          do jo = 1,no
            Haux(1,2,jo,1,io) =  Haux(1,1,io,2,jo)
            Haux(2,2,jo,1,io) = -Haux(2,1,io,2,jo)
          enddo
        enddo
        call cdiag(Haux,Saux,2*nuotot,2*nuotot,2*nuo,eo,psi,
     .             neigwanted,iscf,ierror)

      endif

C Check if we are done ................................................
      if (.not.getD) return

C Find new Fermi energy and occupation weights ........................
      call fermid( 2, 2, 1, (/1.d0/), maxo, neigwanted, eo, 
     .             temp, qtot, qo, ef, Entropy )

*     write(6,'(/,a,/,(10f7.2))') 'diag2g: eo =', eo
*     write(6,'(/,a,/,(10f7.2))') 'diag2g: qo =', qo

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do io = 1,nuotot*2
          qo(io) =  ( stepf( (eo(io)-e2)/t ) -
     .                stepf( (eo(io)-e1)/t ) ) 
        enddo
      endif
      
C New density and energy-density matrices of unit-cell orbitals .......
      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:4) = 0.0d0
      Enew(1:nd,1:4) = 0.0d0

      BNode = 0
      iie = 0
      do ie = 1,nuotot
        if (Node.eq.BNode) then
          iie = iie + 1
        endif
        qe = qo(ie)
        if (abs(qe).gt.occtol) then
          if (Node.eq.BNode) then
            do j = 1,nuotot
              aux(1,j) = psi(1,1,j,iie)
              aux(2,j) = psi(2,1,j,iie)
              aux(3,j) = psi(1,2,j,iie)
              aux(4,j) = psi(2,2,j,iie)
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(aux(1,1),4*no,MPI_double_precision,BNode,
     .      MPI_Comm_World,MPIerror)
#endif
          ee = qo(ie) * eo(ie)
          do io = 1,nuo
            do j = 1,numd(io)
              call LocalToGlobalOrb(io,Node,Nodes,iio)
              ind = listdptr(io) + j
              jo = listd(ind)
              pipj = aux(1,iio) * aux(1,jo) +
     .               aux(2,iio) * aux(2,jo)
              Dnew(ind,1) = Dnew(ind,1) + qe * pipj
              Enew(ind,1) = Enew(ind,1) + ee * pipj
              pipj = aux(3,iio) * aux(3,jo) +
     .               aux(4,iio) * aux(4,jo)
              Dnew(ind,2) = Dnew(ind,2) + qe * pipj
              Enew(ind,2) = Enew(ind,2) + ee * pipj
              pipj = aux(1,iio) * aux(3,jo) +
     .               aux(2,iio) * aux(4,jo)
              Dnew(ind,3) = Dnew(ind,3) + qe * pipj
              Enew(ind,3) = Enew(ind,3) + ee * pipj
              pipj = aux(1,iio) * aux(4,jo) -
     .               aux(2,iio) * aux(3,jo)
              Dnew(ind,4) = Dnew(ind,4) + qe * pipj
              Enew(ind,4) = Enew(ind,4) + ee * pipj
            enddo
          enddo
        endif
        BTest = ie/BlockSize
        if (BTest*BlockSize.eq.ie) then
          BNode = BNode + 1
          if (BNode .gt. Nodes-1) BNode = 0
        endif
      enddo

      end
