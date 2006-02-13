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
      subroutine diag2k( nuo, no, maxuo, maxnh, maxnd, maxo,
     .                   numh, listhptr, listh, numd, listdptr, 
     .                   listd, H, S, getD, qtot, temp, e1, e2,
     .                   xij, indxuo, nk, kpoint, wk,
     .                   eo, qo, Dnew, Enew, ef, Entropy,
     .                   Haux, Saux, psi, Dk, Ek, aux, nuotot,
     .                   occtol, iscf )
C *********************************************************************
C Calculates the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for non-colinear spin with k-sampling and time 
C reversal symmetry.
C Written by J.Soler, August 1998.
C Modified by V.M.Garcia, June 2002.
C Density matrix build and computation of eigenvectors restricted
C by J.D. Gale November 2004.
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer maxuo               : Maximum number of basis orbitals
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd / DM
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
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C real*8  occtol              : Occupancy threshold for DM build
C integer iscf                : SCF cycle number
C *************************** OUTPUT **********************************
C real*8 eo(maxo*4,nk)        : Eigenvalues
C real*8 qo(maxo*4,nk)        : Occupations of eigenstates
C real*8 Dnew(maxnd,4)        : Output Density Matrix
C real*8 Enew(maxnd,4)        : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,2,nuotot,2,nuo) : Aux. space for the hamiltonian matrix
C real*8 Saux(2,2,nuotot,2,nuo) : Aux. space for the overlap matrix
C real*8 psi(2,2,nuotot,2*nuo)  : Aux. space for the eigenvectors
C real*8 aux(5,2*nuotot)        : Extra auxiliary space
C real*8 Dk(2,2,nuotot,2,nuo)   : Aux. space that may be the same as Haux
C real*8 Ek(2,2,nuotot,2,nuo)   : Aux. space that may be the same as Saux
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
      use parallel,      only : Node, Nodes, BlockSize
      use parallelsubs,  only : LocalToGlobalOrb
      use m_fermid,      only : fermid, stepf
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer
     .  maxuo, maxnd, maxnh, maxo, nk, no, nuo, nuotot, iscf

      integer 
     .  indxuo(no), listh(maxnh), numh(nuo), listd(maxnd), numd(nuo),
     .  listhptr(*), listdptr(*)

      real(dp)
     .  Dnew(maxnd,4), occtol,
     .  e1, e2, ef, Enew(maxnd,4), eo(maxo*4,nk), Entropy,
     .  H(maxnh,4), kpoint(3,nk), qo(maxo*4,nk), qtot,
     .  S(maxnh), temp, wk(nk), xij(3,maxnh)
     
      real(dp)
     .  aux(5,2*nuotot), Dk(2,2,nuotot,2,nuo), Ek(2,2,nuotot,2,nuo), 
     .  Haux(2,2,nuotot,2,nuo), psi(2,2,nuotot,nuo*2), 
     .  Saux(2,2,nuotot,2,nuo)

      logical
     .  getD

      external    cdiag

C  Internal variables .............................................
      integer
     .  BNode, BTest, ie, ierror, iie, ik, ind, io, iio,
     .  iuo, j, jo, juo, nd, neigneeded
      real(dp)
     .  ckx, ee, kxij, pipj1, pipj2, qe, skx, t
C  ....................

C Find eigenvalues at every k point ...............................
      do ik = 1,nk

C       Initialize Hamiltonian and overlap matrices in full format
C       Index i is for real/imag parts
C       Indices is and js are for spin components
C       Indices iuo and juo are for orbital components:
C       Haux(i,js,juo,is,iuo) = <js,juo|H|is,iuo>
        Saux = 0.0d0
        Haux = 0.0d0

C       Transfer S,H matrices from sparse format in supercell to
C       full format in unit cell
C       Convention: ispin=1 => H11, ispin=2 => H22, 
C                   ispin=3 => Real(H12), ispin=4 => Imag(H12)
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,ind) +
     .             kpoint(2,ik) * xij(2,ind) +
     .             kpoint(3,ik) * xij(3,ind)
            ckx = cos(kxij)
            skx = sin(kxij)
            Saux(1,1,juo,1,iuo) = Saux(1,1,juo,1,iuo) + S(ind)*ckx
            Saux(2,1,juo,1,iuo) = Saux(2,1,juo,1,iuo) + S(ind)*skx
            Saux(1,2,juo,2,iuo) = Saux(1,2,juo,2,iuo) + S(ind)*ckx
            Saux(2,2,juo,2,iuo) = Saux(2,2,juo,2,iuo) + S(ind)*skx
            Haux(1,1,juo,1,iuo) = Haux(1,1,juo,1,iuo) + H(ind,1)*ckx
            Haux(2,1,juo,1,iuo) = Haux(2,1,juo,1,iuo) + H(ind,1)*skx
            Haux(1,2,juo,2,iuo) = Haux(1,2,juo,2,iuo) + H(ind,2)*ckx
            Haux(2,2,juo,2,iuo) = Haux(2,2,juo,2,iuo) + H(ind,2)*skx
            Haux(1,1,juo,2,iuo) = Haux(1,1,juo,2,iuo) + H(ind,3)*ckx
     .                                                + H(ind,4)*skx
            Haux(2,1,juo,2,iuo) = Haux(2,1,juo,2,iuo) - H(ind,4)*ckx
     .                                                + H(ind,3)*skx
          enddo
        enddo

C Hermiticity
        do iuo = 1,nuo
          do juo = 1,nuo
            Haux(1,2,juo,1,iuo) =   Haux(1,1,iuo,2,juo)
            Haux(2,2,juo,1,iuo) = - Haux(2,1,iuo,2,juo)
          enddo
        enddo

C Find eigenvalues
C Possible memory optimization: equivalence Haux and psi
        call cdiag(Haux,Saux,2*nuotot,2*nuotot,2*nuo,eo(1,ik),psi,
     .             0,iscf,ierror)
        if (ierror.ne.0) then
          call die('Terminating due to failed diagonalisation')
        endif
      enddo

C Check if we are done ................................................
      if (.not.getD) return

C Find new Fermi energy and occupation weights ........................
      call fermid( 2, 4, nk, wk, maxo, nuotot, eo, 
     .             temp, qtot, qo, ef, Entropy )

*     write(6,'(/,a,/,(10f7.2))') 'diag2k: eo =',(eo(ie,1),ie=1,nuo*2)
*     write(6,'(/,a,/,(10f7.2))') 'diag2k: qo =',(qo(ie,1),ie=1,nuo*2)
*     write(6,'(/,a)') 'diag2k: eo ='
*     do ik = 1,nk
*       write(6,'(a,i6,/,(10f7.2))') 'ik=',ik,(eo(ie,ik),ie=1,nuo*2)
*     enddo

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ik = 1,nk
          do io = 1,nuotot*2
            qo(io,ik) = wk(ik) * 
     .           ( stepf( (eo(io,ik)-e2)/t ) -
     .             stepf( (eo(io,ik)-e1)/t ) ) 
          enddo
        enddo
      endif
      
C New density and energy-density matrices of unit-cell orbitals .......
      nd = listdptr(nuo) + numd(nuo)
      Dnew(1:nd,1:4) = 0.0d0
      Enew(1:nd,1:4) = 0.0d0

      do ik = 1,nk

C Find maximum eigenvector that is required for this k point 
        neigneeded = 0
        ie = nuotot
        do while (ie.gt.0.and.neigneeded.eq.0)
          qe = qo(ie,ik)
          if (qe.gt.occtol) neigneeded = ie
          ie = ie - 1
        enddo

C Find eigenvectors 
        Saux = 0.0d0
        Haux = 0.0d0
        do iuo = 1,nuo
          do j = 1,numh(iuo)
            ind = listhptr(iuo) + j
            jo = listh(ind)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,ind) +
     .             kpoint(2,ik) * xij(2,ind) +
     .             kpoint(3,ik) * xij(3,ind)
            ckx = cos(kxij)
            skx = sin(kxij)
            Saux(1,1,juo,1,iuo) = Saux(1,1,juo,1,iuo) + S(ind)*ckx
            Saux(2,1,juo,1,iuo) = Saux(2,1,juo,1,iuo) + S(ind)*skx
            Saux(1,2,juo,2,iuo) = Saux(1,2,juo,2,iuo) + S(ind)*ckx
            Saux(2,2,juo,2,iuo) = Saux(2,2,juo,2,iuo) + S(ind)*skx
            Haux(1,1,juo,1,iuo) = Haux(1,1,juo,1,iuo) + H(ind,1)*ckx
            Haux(2,1,juo,1,iuo) = Haux(2,1,juo,1,iuo) + H(ind,1)*skx
            Haux(1,2,juo,2,iuo) = Haux(1,2,juo,2,iuo) + H(ind,2)*ckx
            Haux(2,2,juo,2,iuo) = Haux(2,2,juo,2,iuo) + H(ind,2)*skx
            Haux(1,1,juo,2,iuo) = Haux(1,1,juo,2,iuo) + H(ind,3)*ckx
     .                                                + H(ind,4)*skx
            Haux(2,1,juo,2,iuo) = Haux(2,1,juo,2,iuo) - H(ind,4)*ckx
     .                                                + H(ind,3)*skx
          enddo
        enddo

C Hermiticity
        do iuo = 1,nuo
          do juo = 1,nuo
            Haux(1,2,juo,1,iuo) =   Haux(1,1,iuo,2,juo)
            Haux(2,2,juo,1,iuo) = - Haux(2,1,iuo,2,juo)
          enddo
        enddo

        call cdiag(Haux,Saux,2*nuotot,2*nuotot,2*nuo,aux,psi,
     .             neigneeded,iscf,ierror)

C Check error flag and take appropriate action
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
C Repeat diagonalisation with increased memory to handle clustering
          Saux = 0.0d0
          Haux = 0.0d0
          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xij(1,ind) +
     .               kpoint(2,ik) * xij(2,ind) +
     .               kpoint(3,ik) * xij(3,ind)
              ckx = cos(kxij)
              skx = sin(kxij)
              Saux(1,1,juo,1,iuo) = Saux(1,1,juo,1,iuo) + S(ind)*ckx
              Saux(2,1,juo,1,iuo) = Saux(2,1,juo,1,iuo) + S(ind)*skx
              Saux(1,2,juo,2,iuo) = Saux(1,2,juo,2,iuo) + S(ind)*ckx
              Saux(2,2,juo,2,iuo) = Saux(2,2,juo,2,iuo) + S(ind)*skx
              Haux(1,1,juo,1,iuo) = Haux(1,1,juo,1,iuo) + H(ind,1)*ckx
              Haux(2,1,juo,1,iuo) = Haux(2,1,juo,1,iuo) + H(ind,1)*skx
              Haux(1,2,juo,2,iuo) = Haux(1,2,juo,2,iuo) + H(ind,2)*ckx
              Haux(2,2,juo,2,iuo) = Haux(2,2,juo,2,iuo) + H(ind,2)*skx
              Haux(1,1,juo,2,iuo) = Haux(1,1,juo,2,iuo) + H(ind,3)*ckx
     .                                                  + H(ind,4)*skx
              Haux(2,1,juo,2,iuo) = Haux(2,1,juo,2,iuo) - H(ind,4)*ckx
     .                                                  + H(ind,3)*skx
            enddo
          enddo
          do iuo = 1,nuo
            do juo = 1,nuo
              Haux(1,2,juo,1,iuo) =  Haux(1,1,iuo,2,juo)
              Haux(2,2,juo,1,iuo) = -Haux(2,1,iuo,2,juo)
            enddo
          enddo
          call cdiag(Haux,Saux,2*nuotot,2*nuotot,2*nuo,aux,psi,
     .               neigneeded,iscf,ierror)
        endif

C Store the products of eigenvectors in matrices Dk and Ek
C WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
        Dk = 0.0d0
        Ek = 0.0d0

        BNode = 0
        iie = 0
        do ie = 1,nuotot
          if (Node.eq.BNode) then
            iie = iie + 1
          endif
          qe = qo(ie,ik)
          if (qe.gt.occtol) then
            if (Node.eq.BNode) then
              do j = 1,nuotot
                aux(1,j) = psi(1,1,j,iie)
                aux(2,j) = psi(2,1,j,iie)
                aux(3,j) = psi(1,2,j,iie)
                aux(4,j) = psi(2,2,j,iie)
              enddo
            endif
#ifdef MPI
            call MPI_Bcast(aux(1,1),4*nuotot,MPI_double_precision,BNode,
     .        MPI_Comm_World,MPIerror)
#endif
            ee = qo(ie,ik) * eo(ie,ik)
            do iuo = 1,nuo
              call LocalToGlobalOrb(iuo,Node,Nodes,iio)
              do juo = 1,nuotot

                pipj1 = aux(1,iio) * aux(1,juo) +
     .                  aux(2,iio) * aux(2,juo)
                pipj2 = aux(1,iio) * aux(2,juo) -
     .                  aux(2,iio) * aux(1,juo)
                Dk(1,1,juo,1,iuo) = Dk(1,1,juo,1,iuo) + qe * pipj1
                Dk(2,1,juo,1,iuo) = Dk(2,1,juo,1,iuo) + qe * pipj2
                Ek(1,1,juo,1,iuo) = Ek(1,1,juo,1,iuo) + ee * pipj1
                Ek(2,1,juo,1,iuo) = Ek(2,1,juo,1,iuo) + ee * pipj2

                pipj1 = aux(3,iio) * aux(3,juo) +
     .                  aux(4,iio) * aux(4,juo)
                pipj2 = aux(3,iio) * aux(4,juo) -
     .                  aux(4,iio) * aux(3,juo)
                Dk(1,2,juo,2,iuo) = Dk(1,2,juo,2,iuo) + qe * pipj1
                Dk(2,2,juo,2,iuo) = Dk(2,2,juo,2,iuo) + qe * pipj2
                Ek(1,2,juo,2,iuo) = Ek(1,2,juo,2,iuo) + ee * pipj1
                Ek(2,2,juo,2,iuo) = Ek(2,2,juo,2,iuo) + ee * pipj2

                pipj1 = aux(1,iio) * aux(3,juo) +
     .                  aux(2,iio) * aux(4,juo)
                pipj2 = aux(1,iio) * aux(4,juo) -
     .                  aux(2,iio) * aux(3,juo)
                Dk(1,1,juo,2,iuo) = Dk(1,1,juo,2,iuo) + qe * pipj1
                Dk(2,1,juo,2,iuo) = Dk(2,1,juo,2,iuo) + qe * pipj2
                Ek(1,1,juo,2,iuo) = Ek(1,1,juo,2,iuo) + ee * pipj1
                Ek(2,1,juo,2,iuo) = Ek(2,1,juo,2,iuo) + ee * pipj2

              enddo
            enddo
          endif
          BTest = ie/BlockSize
          if (BTest*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
        enddo

C Add contribution to density matrices of unit-cell orbitals
        do iuo = 1,nuo
          do j = 1,numd(iuo)
            ind = listdptr(iuo) + j
            jo = listd(ind)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,ind) +
     .             kpoint(2,ik) * xij(2,ind) +
     .             kpoint(3,ik) * xij(3,ind)
            ckx = cos(kxij)
            skx = sin(kxij)

            Dnew(ind,1) = Dnew(ind,1) + Dk(1,1,juo,1,iuo) * ckx
     .                                + Dk(2,1,juo,1,iuo) * skx
            Enew(ind,1) = Enew(ind,1) + Ek(1,1,juo,1,iuo) * ckx
     .                                + Ek(2,1,juo,1,iuo) * skx

            Dnew(ind,2) = Dnew(ind,2) + Dk(1,2,juo,2,iuo) * ckx
     .                                + Dk(2,2,juo,2,iuo) * skx
            Enew(ind,2) = Enew(ind,2) + Ek(1,2,juo,2,iuo) * ckx
     .                                + Ek(2,2,juo,2,iuo) * skx

C Average k and -k solutions because time-reversal sym
            Dnew(ind,3) = Dnew(ind,3) + 0.5d0 * 
     .       ((Dk(1,1,juo,2,iuo) + Dk(1,1,iuo,2,juo)) * ckx
     .       +(Dk(2,1,juo,2,iuo) - Dk(2,1,iuo,2,juo)) * skx)
            Enew(ind,3) = Enew(ind,3) + 0.5d0 * 
     .       ((Ek(1,1,juo,2,iuo) + Ek(1,1,iuo,2,juo)) * ckx
     .       +(Ek(2,1,juo,2,iuo) - Ek(2,1,iuo,2,juo)) * skx)

            Dnew(ind,4) = Dnew(ind,4) + 0.5d0 * 
     .       ((Dk(2,1,juo,2,iuo) + Dk(2,1,iuo,2,juo)) * ckx
     .       -(Dk(1,1,juo,2,iuo) - Dk(1,1,iuo,2,juo)) * skx)
            Enew(ind,4) = Enew(ind,4) + 0.5d0 * 
     .       ((Ek(2,1,juo,2,iuo) + Ek(2,1,iuo,2,juo)) * ckx
     .       -(Ek(1,1,juo,2,iuo) - Ek(1,1,iuo,2,juo)) * skx)

          enddo
        enddo

      enddo

      end
