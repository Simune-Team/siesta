      subroutine diag2k( nuo, no, maxuo, maxno, maxnd, maxo,
     .                   numh, listh, numd, listd, H, S,
     .                   getD, qtot, temp, e1, e2,
     .                   xij, indxuo, nk, kpoint, wk,
     .                   eo, qo, Dnew, Enew, ef,
     .                   Haux, Saux, psi, Dk, Ek, aux,
     .                   nuotot, Node, Nodes)
C *********************************************************************
C Calculates the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for non-colinear spin with k-sampling.
C Writen by J.Soler, August 1998.
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals in unit cell
C integer no                  : Number of basis orbitals in supercell
C integer maxuo               : Maximum number of basis orbitals
C integer maxno               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
C integer numh(no)            : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listh(maxno,no)     : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(no)            : Number of nonzero elements of each row 
C                               of density matrix
C integer listd(maxnd,no)     : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxno,maxuo,4)    : Hamiltonian in sparse form
C real*8  S(maxno,maxuo)      : Overlap in sparse form
C logical getD                : Find occupations and density matrices?
C real*8  qtot                : Number of electrons in unit cell
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
C real*8  xij(3,maxno,maxuo)  : Vectors between orbital centers (sparse)
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
C integer Node                : local processor
C integer Nodes               : number of processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo*4,nk)      : Eigenvalues
C real*8 qo(maxo*4,nk)      : Occupations of eigenstates
C real*8 Dnew(maxnd,maxuo,4): Output Density Matrix
C real*8 Enew(maxnd,maxuo,4): Output Energy-Density Matrix
C real*8 ef                 : Fermi energy
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
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

#ifdef MPI
      integer 
     .  MPIerror
#endif

      integer
     .  maxuo, maxnd, maxno, maxo, nk, no, nuo, nuotot, Node, Nodes

      integer 
     .  indxuo(no), listh(maxno,no), numh(no),
     .  listd(maxnd,no), numd(no)

      double precision
     .  Dnew(maxnd,maxuo,4),
     .  e1, e2, ef, Enew(maxnd,maxuo,4), eo(maxo*4,nk),
     .  H(maxno,maxuo,4), kpoint(3,nk), qo(maxo*4,nk), qtot,
     .  S(maxno,maxuo), stepf, temp, wk(nk), xij(3,maxno,maxuo)
     
      double precision
     .  aux(5,2*nuotot), Dk(2,2,nuotot,2,nuo), Ek(2,2,nuotot,2,nuo), 
     .  Haux(2,2,nuotot,2,nuo), psi(2,2,nuotot,nuo*2), 
     .  Saux(2,2,nuotot,2,nuo)

      logical
     .  getD

      external
     .  cdiag, fermid, stepf

C  Internal variables .............................................
      integer
     .  i, ie, ik, io, iio, is, ispin, iuo, j, jo, js, juo, 
     .  BNode, iie, BTest, ierror
      double precision
     .  ckx, ee, kxij, pipj1, pipj2, qe, skx, t
C  ....................

C Find eigenvalues at every k point ...............................
      do ik = 1,nk

C       Initialize Hamiltonian and overlap matrices in full format
C       Index i is for real/imag parts
C       Indices is and js are for spin components
C       Indices iuo and juo are for orbital components:
C       Haux(i,js,juo,is,iuo) = <js,juo|H|is,iuo>
   10   do iuo = 1,nuo
          do is = 1,2
            do juo = 1,nuotot
              do js = 1,2
                do i = 1,2
                  Saux(i,js,juo,is,iuo) = 0.d0
                  Haux(i,js,juo,is,iuo) = 0.d0
                enddo
              enddo
            enddo
          enddo
        enddo

C       Transfer S,H matrices from sparse format in supercell to
C       full format in unit cell
C       Convention: ispin=1 => H11, ispin=2 => H22, 
C                   ispin=3 => Real(H12), ispin=4 => Imag(H12)
        do io = 1,nuo
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do j = 1,numh(iio)
            jo = listh(j,iio)
            call GlobalToLocalOrb(indxuo(iio),Node,Nodes,iuo)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,j,io) +
     .             kpoint(2,ik) * xij(2,j,io) +
     .             kpoint(3,ik) * xij(3,j,io)
            ckx = cos(kxij)
            skx = sin(kxij)
            Saux(1,1,juo,1,iuo) = Saux(1,1,juo,1,iuo) + S(j,io)*ckx
            Saux(2,1,juo,1,iuo) = Saux(2,1,juo,1,iuo) + S(j,io)*skx
            Saux(1,2,juo,2,iuo) = Saux(1,2,juo,2,iuo) + S(j,io)*ckx
            Saux(2,2,juo,2,iuo) = Saux(2,2,juo,2,iuo) + S(j,io)*skx
            Haux(1,1,juo,1,iuo) = Haux(1,1,juo,1,iuo) + H(j,io,1)*ckx
            Haux(2,1,juo,1,iuo) = Haux(2,1,juo,1,iuo) + H(j,io,1)*skx
            Haux(1,2,juo,2,iuo) = Haux(1,2,juo,2,iuo) + H(j,io,2)*ckx
            Haux(2,2,juo,2,iuo) = Haux(2,2,juo,2,iuo) + H(j,io,2)*skx
            Haux(1,1,juo,2,iuo) = Haux(1,1,juo,2,iuo) + H(j,io,3)*ckx
     .                                                + H(j,io,4)*skx
            Haux(2,1,juo,2,iuo) = Haux(2,1,juo,2,iuo) - H(j,io,4)*ckx
     .                                                + H(j,io,3)*skx
            Haux(1,2,juo,1,iuo) = Haux(1,2,juo,1,iuo) + H(j,io,3)*ckx
     .                                                - H(j,io,4)*skx
            Haux(2,2,juo,1,iuo) = Haux(2,2,juo,1,iuo) + H(j,io,4)*ckx
     .                                                + H(j,io,3)*skx
          enddo
        enddo

C       Solve the eigenvalue problem 
C       Possible CPU-time optimization: only eigenvalues, not eigvect
C       Possible memory optimization: equivalence Haux and psi
        call cdiag( Haux, 2*nuotot, Saux, 2*nuotot, 2*nuo,
     .              eo(1,ik), psi, 2*nuotot, aux, ierror )
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
          goto 10
        endif
      enddo
C ....................

C Check if we are done ................................................
      if (.not.getD) return
C ....................

C Find new Fermi energy and occupation weights ........................
      call fermid( 2, 4, nk, wk, maxo, nuotot, eo, 
     .             temp, qtot, qo, ef )
C ....................

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
     .             stepf( (eo(io,ik)-e1)/t ) ) / 2
          enddo
        enddo
      endif
C ....................
      
c New density and energy-density matrices of unit-cell orbitals .......
      do ispin = 1,4
        do io = 1,nuo
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do j = 1,numd(iio)
            Dnew(j,io,ispin) = 0.d0
            Enew(j,io,ispin) = 0.d0
          enddo
        enddo
      enddo

      do ik = 1,nk

C       Find eigenvectors again (were stored only for one k point)
   20   do iuo = 1,nuo
          do is = 1,2
            do juo = 1,nuotot
              do js = 1,2
                do i = 1,2
                  Saux(i,js,juo,is,iuo) = 0.d0
                  Haux(i,js,juo,is,iuo) = 0.d0
                enddo
              enddo
            enddo
          enddo
        enddo
        do io = 1,nuo
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do j = 1,numh(iio)
            jo = listh(j,iio)
            call GlobalToLocalOrb(indxuo(iio),Node,Nodes,iuo)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,j,io) +
     .             kpoint(2,ik) * xij(2,j,io) +
     .             kpoint(3,ik) * xij(3,j,io)
            ckx = cos(kxij)
            skx = sin(kxij)
            Saux(1,1,juo,1,iuo) = Saux(1,1,juo,1,iuo) + S(j,io)*ckx
            Saux(2,1,juo,1,iuo) = Saux(2,1,juo,1,iuo) + S(j,io)*skx
            Saux(1,2,juo,2,iuo) = Saux(1,2,juo,2,iuo) + S(j,io)*ckx
            Saux(2,2,juo,2,iuo) = Saux(2,2,juo,2,iuo) + S(j,io)*skx
            Haux(1,1,juo,1,iuo) = Haux(1,1,juo,1,iuo) + H(j,io,1)*ckx
            Haux(2,1,juo,1,iuo) = Haux(2,1,juo,1,iuo) + H(j,io,1)*skx
            Haux(1,2,juo,2,iuo) = Haux(1,2,juo,2,iuo) + H(j,io,2)*ckx
            Haux(2,2,juo,2,iuo) = Haux(2,2,juo,2,iuo) + H(j,io,2)*skx
            Haux(1,1,juo,2,iuo) = Haux(1,1,juo,2,iuo) + H(j,io,3)*ckx
     .                                                + H(j,io,4)*skx
            Haux(2,1,juo,2,iuo) = Haux(2,1,juo,2,iuo) - H(j,io,4)*ckx
     .                                                + H(j,io,3)*skx
            Haux(1,2,juo,1,iuo) = Haux(1,2,juo,1,iuo) + H(j,io,3)*ckx
     .                                                - H(j,io,4)*skx
            Haux(2,2,juo,1,iuo) = Haux(2,2,juo,1,iuo) + H(j,io,4)*ckx
     .                                                + H(j,io,3)*skx
          enddo
        enddo
        call cdiag( Haux, 2*nuotot, Saux, 2*nuotot, 2*nuo,
     .              eo(1,ik), psi, 2*nuotot, aux, ierror )
        if (ierror.gt.0) then
          call die('Terminating due to failed diagonalisation')
        elseif (ierror.lt.0) then
          goto 20
        endif

C       Store the products of eigenvectors in matrices Dk and Ek
C       WARNING: Dk and Ek may be EQUIVALENCE'd to Haux and Saux
        do iuo = 1,nuo
          do is = 1,2
            do juo = 1,nuotot
              do js = 1,2
                do i = 1,2
                  Dk(i,js,juo,is,iuo) = 0.d0
                  Ek(i,js,juo,is,iuo) = 0.d0
                enddo
              enddo
            enddo
          enddo
        enddo
        BNode = 0
        iie = 0
        do ie = 1,nuotot*2
          if (Node.eq.BNode) then
            iie = iie + 1
            do j = 1,nuotot
              aux(1,j) = psi(1,1,j,iie)
              aux(2,j) = psi(2,1,j,iie)
              aux(3,j) = psi(1,2,j,iie)
              aux(4,j) = psi(2,2,j,iie)
            enddo
          endif
#ifdef MPI
          call MPI_Bcast(aux(1,1),4*nuotot,DAT_double,BNode,
     .      MPI_Comm_World,MPIerror)
#endif
          qe = qo(ie,ik)
          ee = qo(ie,ik) * eo(ie,ik)
          do iuo = 1,nuo
            call LocalToGlobalOrb(iuo,Node,Nodes,iio)
            do juo = 1,nuotot

              pipj1 = aux(1,iio) * aux(1,juo) +
     .                aux(2,iio) * aux(2,juo)
              pipj2 = aux(1,iio) * aux(2,juo) -
     .                aux(2,iio) * aux(1,juo)
              Dk(1,1,juo,1,iuo) = Dk(1,1,juo,1,iuo) + qe * pipj1
              Dk(2,1,juo,1,iuo) = Dk(2,1,juo,1,iuo) + qe * pipj2
              Ek(1,1,juo,1,iuo) = Ek(1,1,juo,1,iuo) + ee * pipj1
              Ek(2,1,juo,1,iuo) = Ek(2,1,juo,1,iuo) + ee * pipj2

              pipj1 = aux(3,iio) * aux(3,juo) +
     .                aux(4,iio) * aux(4,juo)
              pipj2 = aux(3,iio) * aux(4,juo) -
     .                aux(4,iio) * aux(3,juo)
              Dk(1,2,juo,2,iuo) = Dk(1,2,juo,2,iuo) + qe * pipj1
              Dk(2,2,juo,2,iuo) = Dk(2,2,juo,2,iuo) + qe * pipj2
              Ek(1,2,juo,2,iuo) = Ek(1,2,juo,2,iuo) + ee * pipj1
              Ek(2,2,juo,2,iuo) = Ek(2,2,juo,2,iuo) + ee * pipj2

              pipj1 = aux(1,iio) * aux(3,juo) +
     .                aux(2,iio) * aux(4,juo)
              pipj2 = aux(1,iio) * aux(4,juo) -
     .                aux(2,iio) * aux(3,juo)
              Dk(1,1,juo,2,iuo) = Dk(1,1,juo,2,iuo) + qe * pipj1
              Dk(2,1,juo,2,iuo) = Dk(2,1,juo,2,iuo) + qe * pipj2
              Ek(1,1,juo,2,iuo) = Ek(1,1,juo,2,iuo) + ee * pipj1
              Ek(2,1,juo,2,iuo) = Ek(2,1,juo,2,iuo) + ee * pipj2

*             pipj1 = psi(1,2,iio,ie) * psi(1,1,juo,ie) +
*    .                psi(2,2,iio,ie) * psi(2,1,juo,ie)
*             pipj2 = psi(1,2,iio,ie) * psi(2,1,juo,ie) -
*    .                psi(2,2,iio,ie) * psi(1,1,juo,ie)
*             Dk(1,2,juo,1,iuo) = Dk(1,2,juo,1,iuo) + qe * pipj1
*             Dk(2,2,juo,1,iuo) = Dk(2,2,juo,1,iuo) + qe * pipj2
*             Ek(1,2,juo,1,iuo) = Ek(1,2,juo,1,iuo) + ee * pipj1
*             Ek(2,2,juo,1,iuo) = Ek(2,2,juo,1,iuo) + ee * pipj2
            enddo
          enddo
          BTest = ie/BlockSize
          if (BTest*BlockSize.eq.ie) then
            BNode = BNode + 1
            if (BNode .gt. Nodes-1) BNode = 0
          endif
        enddo

C       Add contribution to density matrices of unit-cell orbitals
        do io = 1,nuo
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          call GlobalToLocalOrb(indxuo(iio),Node,Nodes,iuo)
          do j = 1,numd(iio)
            jo = listd(j,iio)
            juo = indxuo(jo)
            kxij = kpoint(1,ik) * xij(1,j,io) +
     .             kpoint(2,ik) * xij(2,j,io) +
     .             kpoint(3,ik) * xij(3,j,io)
            ckx = cos(kxij)
            skx = sin(kxij)
            Dnew(j,io,1) = Dnew(j,io,1) + Dk(1,1,juo,1,iuo) * ckx
     .                                  + Dk(2,1,juo,1,iuo) * skx
            Enew(j,io,1) = Enew(j,io,1) + Ek(1,1,juo,1,iuo) * ckx
     .                                  + Ek(2,1,juo,1,iuo) * skx
            Dnew(j,io,2) = Dnew(j,io,2) + Dk(1,2,juo,2,iuo) * ckx
     .                                  + Dk(2,2,juo,2,iuo) * skx
            Enew(j,io,2) = Enew(j,io,2) + Ek(1,2,juo,2,iuo) * ckx
     .                                  + Ek(2,2,juo,2,iuo) * skx
            Dnew(j,io,3) = Dnew(j,io,3) + Dk(1,1,juo,2,iuo) * ckx
     .                                  + Dk(2,1,juo,2,iuo) * skx
            Enew(j,io,3) = Enew(j,io,3) + Ek(1,1,juo,2,iuo) * ckx
     .                                  + Ek(2,1,juo,2,iuo) * skx
            Dnew(j,io,4) = Dnew(j,io,4) + Dk(2,1,juo,2,iuo) * ckx
     .                                  - Dk(1,1,juo,2,iuo) * skx
            Enew(j,io,4) = Enew(j,io,4) + Ek(2,1,juo,2,iuo) * ckx
     .                                  - Ek(1,1,juo,2,iuo) * skx
          enddo
        enddo

      enddo
C ....................

      end


