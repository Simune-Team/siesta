      subroutine diag2g( nuo, no, maxuo, maxno, maxnd, maxo,
     .                   numh, listh, numd, listd, H, S,
     .                   getD, qtot, temp, e1, e2,
     .                   eo, qo, Dnew, Enew, ef,
     .                   Haux, Saux, psi, aux, nuotot,
     .                   Node, Nodes)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices.
C This version is for non-colinear spin at gamma point.
C Writen by J.Soler, May and August 1998.
C **************************** INPUT **********************************
C integer nuo                 : Number of basis orbitals on local node
C integer no                  : Number of basis orbitals
C integer maxuo               : Maximum number of basis  orbitals
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
C                               ofdensity matrix
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
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C integer Node                : local processor
C integer Nodes               : number of processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo*2)          : Eigenvalues
C real*8 qo(maxo*2)          : Occupations of eigenstates
C real*8 Dnew(maxnd,maxuo,4) : Output Density Matrix
C real*8 Enew(maxnd,maxuo,4) : Output Energy-Density Matrix
C real*8 ef                     : Fermi energy
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
     .  maxo, maxuo, maxnd, maxno, no, nuo, nuotot, Node, Nodes

      integer 
     .  listh(maxno,no), numh(no),
     .  listd(maxnd,no), numd(no)

      double precision
     .  Dnew(maxnd,maxuo,4),
     .  e1, e2, ef, Enew(maxnd,maxuo,4), eo(maxo*2),
     .  H(maxno,maxuo,4), qo(maxo*2), qtot,
     .  S(maxno,maxo), stepf, temp
     
      double precision
     .  aux(5,2*nuotot), Haux(2,2,nuotot,2,nuo),
     .  psi(2,2,nuotot,2*no), Saux(2,2,nuotot,2,nuo)

      logical
     .  getD

      external
     .  cdiag, fermid, stepf

C  Internal variables .............................................
      integer           i, ie, io, iio, is, ispin, j, jo, js, 
     .                  BNode, iie, BTest, ierror
      double precision  ee, pipj, qe, t
C  ....................

C Initialize Hamiltonian and overlap matrices in full format .....
C Index i is for real/imag parts
C Indices is and js are for spin components
C Indices iuo and juo are for orbital components:
C Haux(i,js,juo,is,iuo) = <js,juo|H|is,iuo>
   10 do io = 1,nuo
        do is = 1,2
          do jo = 1,nuotot
            do js = 1,2
              do i = 1,2
                Saux(i,js,jo,is,io) = 0.d0
                Haux(i,js,jo,is,io) = 0.d0
              enddo
            enddo
          enddo
        enddo
      enddo
C  ....................

C Copy S,H matrices from sparse to full format ....................
C Convention: ispin=1 => H11, ispin=2 => H22, ispin=3 => Real(H12),
C             ispin=4 => Imag(H12)
      do io = 1,nuo
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        do j = 1,numh(iio)
          jo = listh(j,iio)
          Saux(1,1,jo,1,io) =  S(j,io)
          Saux(1,2,jo,2,io) =  S(j,io)
          Haux(1,1,jo,1,io) =  H(j,io,1)
          Haux(1,2,jo,2,io) =  H(j,io,2)
          Haux(1,1,jo,2,io) =  H(j,io,3)
          Haux(1,2,jo,1,io) =  H(j,io,3)
          Haux(2,1,jo,2,io) = -H(j,io,4)
          Haux(2,2,jo,1,io) =  H(j,io,4)
        enddo
      enddo
C ....................

C Solve the eigenvalue problem .......................................
      call cdiag( Haux, 2*nuotot, Saux, 2*nuotot, 2*nuo,
     .            eo, psi, 2*nuotot, aux, ierror )
      if (ierror.gt.0) then
        call die('Terminating due to failed diagonalisation')
      elseif (ierror.lt.0) then
        goto 10
      endif
C ....................

C Check if we are done ................................................
      if (.not.getD) return
C ....................

C Find new Fermi energy and occupation weights ........................
      call fermid( 2, 2, 1, 1.d0, maxo, no, eo, 
     .             temp, qtot, qo, ef )
C ....................

*     write(6,'(/,a,/,(10f7.2))') 'diag2g: eo =', eo
*     write(6,'(/,a,/,(10f7.2))') 'diag2g: qo =', qo

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do io = 1,nuotot*2
          qo(io) =  ( stepf( (eo(io)-e2)/t ) -
     .                stepf( (eo(io)-e1)/t ) ) / 2
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

      BNode = 0
      iie = 0
      do ie = 1,no*2
        if (Node.eq.BNode) then
          iie = iie + 1
          do j = 1,no
            aux(1,j) = psi(1,1,j,iie)
            aux(2,j) = psi(2,1,j,iie)
            aux(3,j) = psi(1,2,j,iie)
            aux(4,j) = psi(2,2,j,iie)
          enddo
        endif
#ifdef MPI
        call MPI_Bcast(aux(1,1),4*no,DAT_double,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
        qe = qo(ie)
        ee = qo(ie) * eo(ie)
        do io = 1,nuo
          call LocalToGlobalOrb(io,Node,Nodes,iio)
          do j = 1,numd(iio)
            jo = listd(j,iio)
            pipj = aux(1,iio) * aux(1,jo) +
     .             aux(2,iio) * aux(2,jo)
            Dnew(j,io,1) = Dnew(j,io,1) + qe * pipj
            Enew(j,io,1) = Enew(j,io,1) + ee * pipj
            pipj = aux(3,iio) * aux(3,jo) +
     .             aux(4,iio) * aux(4,jo)
            Dnew(j,io,2) = Dnew(j,io,2) + qe * pipj
            Enew(j,io,2) = Enew(j,io,2) + ee * pipj
            pipj = aux(1,iio) * aux(3,jo) +
     .             aux(2,iio) * aux(4,jo)
            Dnew(j,io,3) = Dnew(j,io,3) + qe * pipj
            Enew(j,io,3) = Enew(j,io,3) + ee * pipj
            pipj = aux(1,iio) * aux(4,jo) -
     .             aux(2,iio) * aux(3,jo)
            Dnew(j,io,4) = Dnew(j,io,4) + qe * pipj
            Enew(j,io,4) = Enew(j,io,4) + ee * pipj
          enddo
        enddo
        BTest = ie/BlockSize
        if (BTest*BlockSize.eq.ie) then
          BNode = BNode + 1
          if (BNode .gt. Nodes-1) BNode = 0
        endif
      enddo
C ....................

      end
