c
c     Put dummy subroutine here so that compilers do not complain
c     about empty file in serial compilation.
c
      subroutine dummy_sub
      end
c
#ifdef MPI
      subroutine diagkp( nspin, nuo, no, maxspn, maxuo, maxno, maxnd, 
     .                  maxo, numh, listh, numd, listd, H, S,
     .                  getD, fixspin, qtot, qs, temp, e1, e2,
     .                  xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs,
     .                  Haux, Saux, psi, Dk, Ek, aux, nuotot,
     .                  Node, Nodes)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization). K-sampling version. 
C Created from diagk, written by J.Soler.
C Uses parallelisation over K points instead of parallelisation 
C within them.
C Parallel modifications by J.Gale, November 1999.
C **************************** INPUT **********************************
C integer nspin               : Number of spin components (1 or 2)
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxspn              : Second dimension of eo and qo
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxo                : First dimension of eo and qo
C integer maxno               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer numh(no)            : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listh(maxno,no)     : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(no)            : Number of nonzero elements of each row 
C                               ofdensity matrix
C integer listd(maxnd,no)     : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxno,maxuo,nspin): Hamiltonian in sparse form
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
C real*8 eo(maxo,maxspn,nk)      : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,maxspn,nk)      : Occupations of eigenstates
C real*8 Dnew(maxno,maxuo,nspin) : Output Density Matrix
C real*8 Enew(maxno,maxuo,nspin) : Output Energy-Density Matrix
C real*8 ef                      : Fermi energy
C *************************** AUXILIARY *******************************
C real*8 Haux(2,nuotot,nuo) : Auxiliary space for the hamiltonian matrix
C real*8 Saux(2,nuotot,nuo) : Auxiliary space for the overlap matrix
C real*8 psi(2,nuotot,nuo)  : Auxiliary space for the eigenvectors
C real*8 aux(2,nuotot*5)    : Extra auxiliary space
C real*8 Dk(2,nuotot,nuo)   : Aux. space that may be the same as Haux
C real*8 Ek(2,nuotot,nuo)   : Aux. space that may be the same as Saux
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
      use mpi

      implicit          none

      integer 
     .  MPIerror

      integer           maxnd, maxno, maxspn, maxuo, maxo, nk, no,
     .                  nspin, nuo, nuotot, Node, Nodes
      integer           indxuo(no), listh(maxno,no), numh(no),
     .                  listd(maxnd,no), numd(no)
      double precision  Dnew(maxnd,maxuo,nspin),
     .                  e1, e2, ef, efs(nspin), Enew(maxnd,maxuo,nspin),
     .                  eo(maxo,maxspn,nk), H(maxno,maxuo,nspin),
     .                  kpoint(3,nk), qo(maxo,maxspn,nk), qtot, 
     .                  S(maxno,maxuo), stepf, temp, wk(nk),
     .                  xij(3,maxno,maxuo), qs(nspin)
      double precision  Dk(2,nuotot,nuotot), Ek(2,nuotot,nuotot),
     .                  Haux(2,nuotot,nuotot), Saux(2,nuotot,nuotot),
     .                  psi(2,nuotot,nuotot), aux(2,nuotot*5)
      logical           getD, fixspin
      external          cdiag, fermid, fermispin, stepf, memory

C  Internal variables .............................................
      integer
     .  ie, ik, io, iio, ispin, iuo, j, jo, juo, BNode, is,
     .  nuoloc, ierror
      double precision
     .  ckxij, ee, kxij, pipj1, pipj2, qe, skxij, t
      double precision, dimension(:,:), allocatable ::
     .  Snew
      double precision, dimension(:,:,:), allocatable ::
     .  Hnew, Dnewloc, Enewloc, xijloc
C  ....................

C Create new distribution of H and S
      nuoloc = nuotot
      allocate(Snew(maxno,nuotot))
      call memory('A','D',maxno*nuotot,'diagkp')
      allocate(Hnew(maxno,nuotot,nspin))
      call memory('A','D',maxno*nuotot*nspin,'diagkp')
      allocate(xijloc(3,maxno,nuotot))
      call memory('A','D',3*maxno*nuotot,'diagkp')
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do is = 1,nspin
            do jo = 1,numh(io)
              Hnew(jo,io,is) = H(jo,iio,is)
            enddo
          enddo
          do jo = 1,numh(io)
            Snew(jo,io) = S(jo,iio)
          enddo
          do jo = 1,numh(io)
            xijloc(1,jo,io) = xij(1,jo,iio)
            xijloc(2,jo,io) = xij(2,jo,iio)
            xijloc(3,jo,io) = xij(3,jo,iio)
          enddo
        endif
        do is = 1,nspin
          call MPI_Bcast(Hnew(1,io,is),numh(io),DAT_double,BNode,
     .      MPI_Comm_World,MPIerror)
        enddo
        call MPI_Bcast(Snew(1,io),numh(io),DAT_double,BNode,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(xijloc(1,1,io),3*numh(io),DAT_double,BNode,
     .    MPI_Comm_World,MPIerror)
      enddo

C Solve eigenvalue problem .........................................
      do ik = 1+Node,nk,Nodes
        do ispin = 1,nspin
          do iuo = 1,nuoloc
            do juo = 1,nuotot
              Saux(1,juo,iuo) = 0.0d0
              Saux(2,juo,iuo) = 0.0d0
              Haux(1,juo,iuo) = 0.0d0
              Haux(2,juo,iuo) = 0.0d0
            enddo
          enddo
          do io = 1,nuoloc
            do j = 1,numh(io)
              jo = listh(j,io)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,j,io) +
     .               kpoint(2,ik) * xijloc(2,j,io) +
     .               kpoint(3,ik) * xijloc(3,j,io)
              ckxij = cos(kxij)
              skxij = sin(kxij)
C Note : sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+Snew(j,io)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-Snew(j,io)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+Hnew(j,io,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-Hnew(j,io,ispin)*skxij
            enddo
          enddo
          call cdiag( Haux, nuotot, Saux, nuotot, nuoloc,
     .                eo(1,ispin,ik), psi, nuotot, aux, ierror )
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif
        enddo
      enddo
C ....................

C Globalise eigenvalues
      BNode = -1
      do ik = 1,nk
        BNode = BNode + 1
        if (BNode.eq.Nodes) BNode = 0
        call MPI_Bcast(eo(1,1,ik),nuotot*nspin,DAT_double,BNode,
     .    MPI_Comm_World,MPIerror)
      enddo

C Check if we are done ................................................
      if (.not.getD) goto 999
C ....................

C Allocate local copies of Dk and Ek for building matrices
      allocate(Dnewloc(maxno,nuoloc,2))
      call memory('A','D',2*maxno*nuoloc,'diagkp')
      allocate(Enewloc(maxno,nuoloc,2))
      call memory('A','D',2*maxno*nuoloc,'diagkp')

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, 1, 1.d0, maxuo, no, eo,
     .               temp, qs, qo, efs )
      else
        call fermid( nspin, maxspn, nk, wk, maxo, nuotot, eo, 
     .             temp, qtot, qo, ef )
      endif
C ....................

C Find weights for local density of states ............................
      if (e1 .lt. e2) then
*       e1 = e1 - ef
*       e2 = e2 - ef
        t = max( temp, 1.d-6 )
        do ik = 1,nk
          do ispin = 1,nspin
            do io = 1,nuotot
              qo(io,ispin,ik) = wk(ik) * 
     .             ( stepf( (eo(io,ispin,ik)-e2)/t ) -
     .               stepf( (eo(io,ispin,ik)-e1)/t ) ) / nspin
            enddo
          enddo
        enddo
      endif
C ....................
      
c New density and energy-density matrices of unit-cell orbitals .......

      do ispin = 1,nspin
        do io = 1,nuoloc
          do j = 1,numd(io)
            Dnewloc(j,io,ispin) = 0.d0
            Enewloc(j,io,ispin) = 0.d0
          enddo
        enddo
      enddo

      do ik = 1+Node,nk,Nodes
        do ispin = 1,nspin

C         Find eigenvectors again (were stored only for one k point)
          do iuo = 1,nuoloc
            do juo = 1,nuotot
              Saux(1,juo,iuo) = 0.d0
              Saux(2,juo,iuo) = 0.d0
              Haux(1,juo,iuo) = 0.d0
              Haux(2,juo,iuo) = 0.d0
            enddo
          enddo
          do io = 1,nuoloc
            do j = 1,numh(io)
              jo = listh(j,io)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,j,io) +
     .               kpoint(2,ik) * xijloc(2,j,io) +
     .               kpoint(3,ik) * xijloc(3,j,io)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+Snew(j,io)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-Snew(j,io)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+Hnew(j,io,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-Hnew(j,io,ispin)*skxij
            enddo
          enddo
          call cdiag( Haux, nuotot, Saux, nuotot, nuoloc,
     .                eo(1,ispin,ik), psi, nuotot, aux, ierror )
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif

C         Add contribution to density matrices of unit-cell orbitals
          do iuo = 1,nuoloc
            do juo = 1,nuotot
              Dk(1,juo,iuo) = 0.d0
              DK(2,juo,iuo) = 0.d0
              Ek(1,juo,iuo) = 0.d0
              Ek(2,juo,iuo) = 0.d0
            enddo
          enddo
C Global operation to form new density matrix
          do ie = 1,nuotot
            do j = 1,nuotot
              aux(1,j) = psi(1,j,ie)
              aux(2,j) = psi(2,j,ie)
            enddo
            qe = qo(ie,ispin,ik)
            ee = qo(ie,ispin,ik) * eo(ie,ispin,ik)
            do iuo = 1,nuoloc
              do juo = 1,nuotot
                pipj1 = aux(1,iuo) * aux(1,juo) +
     .                  aux(2,iuo) * aux(2,juo)
                pipj2 = aux(1,iuo) * aux(2,juo) -
     .                  aux(2,iuo) * aux(1,juo)
                Dk(1,juo,iuo) = Dk(1,juo,iuo) + qe * pipj1
                Dk(2,juo,iuo) = Dk(2,juo,iuo) + qe * pipj2
                Ek(1,juo,iuo) = Ek(1,juo,iuo) + ee * pipj1
                Ek(2,juo,iuo) = Ek(2,juo,iuo) + ee * pipj2
              enddo
            enddo
          enddo
          do io = 1,nuoloc
            iuo = indxuo(io)
            do j = 1,numd(io)
              jo = listd(j,io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,j,io) +
     .               kpoint(2,ik) * xijloc(2,j,io) +
     .               kpoint(3,ik) * xijloc(3,j,io)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Dnewloc(j,io,ispin)=Dnewloc(j,io,ispin)+ 
     .          Dk(1,juo,iuo)*ckxij - Dk(2,juo,iuo)*skxij
              Enewloc(j,io,ispin)=Enewloc(j,io,ispin)+ 
     .          Ek(1,juo,iuo)*ckxij - Ek(2,juo,iuo)*skxij
            enddo
          enddo

        enddo
      enddo

C Globalise Dnew/Enew
      do io = 1,nuoloc
        call WhichNodeOrb(io,Nodes,BNode)
        call GlobalToLocalOrb(io,BNode,Nodes,iio)
        do ispin = 1,nspin
          call MPI_Reduce(Dnewloc(1,io,ispin),Dnew(1,iio,ispin),
     .      numd(io),DAT_double,MPI_sum,BNode,MPI_Comm_World,MPIerror)
          call MPI_Reduce(Enewloc(1,io,ispin),Enew(1,iio,ispin),
     .      numd(io),DAT_double,MPI_sum,BNode,MPI_Comm_World,MPIerror)
        enddo
      enddo
C ....................

C Deallocate local memory needed only if getD is true
      call memory('D','D',size(Dnewloc),'diagkp')
      deallocate(Dnewloc)
      call memory('D','D',size(Enewloc),'diagkp')
      deallocate(Enewloc)

C Exit point 
  999 continue

C Free local memory
      call memory('D','D',size(xijloc),'diagkp')
      deallocate(xijloc)
      call memory('D','D',size(Snew),'diagkp')
      deallocate(Snew)
      call memory('D','D',size(Hnew),'diagkp')
      deallocate(Hnew)

      end
#endif
