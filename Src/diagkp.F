c
c     Put dummy subroutine here so that compilers do not complain
c     about empty file in serial compilation.
c
      subroutine dummy_sub
      end
c
#ifdef MPI
      subroutine diagkp( nspin, nuo, no, maxspn, maxuo, maxnh, maxnd, 
     .                  maxo, numh, listhptr, listh, numd, listdptr, 
     .                  listd, H, S, getD, fixspin, qtot, qs, temp, 
     .                  e1, e2, xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs, Entropy,
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
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd / DM
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix locally
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix locally
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row 
C                               of density matrix locally
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix locally
C integer listd(maxnd)        : Nonzero density-matrix element column 
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
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
C integer Node                : local processor
C integer Nodes               : number of processors
C *************************** OUTPUT **********************************
C real*8 eo(maxo,maxspn,nk)   : Eigenvalues
C ******************** OUTPUT (only if getD=.true.) *******************
C real*8 qo(maxo,maxspn,nk)   : Occupations of eigenstates
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix
C real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 Entropy              : Electronic entropy
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
C The auxiliary arrays are now no longer symmetric and so the order
C of referencing has been changed in several places to reflect this.
C Note : It is assumed in a couple of places that the sparsity of
C H/S and Dscf are the same (which is the case in the code at present)
C and therefore the only one set of pointer arrays are globalised.
C *********************************************************************
C
C  Modules
C
      use precision
      use sys
      use parallel
      use mpi_siesta

      implicit          none

      integer 
     .  MPIerror

      integer           maxnd, maxnh, maxspn, maxuo, maxo, nk, no,
     .                  nspin, nuo, nuotot, Node, Nodes
      integer           indxuo(no), listh(maxnh), numh(nuo),
     .                  listd(maxnd), numd(nuo), listhptr(nuo),
     .                  listdptr(nuo)
      double precision  Dnew(maxnd,nspin),
     .                  e1, e2, ef, efs(nspin), Enew(maxnd,nspin),
     .                  Entropy, eo(maxo,maxspn,nk), H(maxnh,nspin),
     .                  kpoint(3,nk), qo(maxo,maxspn,nk), qtot, 
     .                  S(maxnh), stepf, temp, wk(nk),
     .                  xij(3,maxnh), qs(nspin)
      double precision  Dk(2,nuotot,nuotot), Ek(2,nuotot,nuotot),
     .                  Haux(2,nuotot,nuotot), Saux(2,nuotot,nuotot),
     .                  psi(2,nuotot,nuotot), aux(2,nuotot*5)
      logical           getD, fixspin
      external          cdiag, fermid, fermispin, stepf, memory

C  Internal variables .............................................
      integer
     .  BNode, ie, ierror, ik, io, iio, ind, is, ispin, iuo,
     .  j, jo, juo, maxnhg, nuog 
      integer, dimension(:), allocatable, save ::
     .  numhg, listhptrg, listhg
      double precision
     .  ckxij, ee, kxij, pipj1, pipj2, qe, skxij, t
      double precision, dimension(:), allocatable, save ::
     .  Snew
      double precision, dimension(:,:), allocatable, save ::
     .  Hnew, Dnewloc, Enewloc, xijloc
C  ....................

C Globalise list arrays - assumes listh and listd are the same

C Allocate local memory for global list arrays
      allocate(numhg(nuotot))
      call memory('A','I',nuotot,'diagkp')
      allocate(listhptrg(nuotot))
      call memory('A','I',nuotot,'diagkp')

C Globalise numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
        endif
        call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
      enddo

C Build global listhptr
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

C Globalse listh
      maxnhg = listhptrg(nuotot)+numhg(nuotot)
      allocate(listhg(maxnhg))
      call memory('A','I',maxnhg,'diagkp')
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = 
     .        listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
          enddo
        endif
        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,
     .    BNode,MPI_Comm_World,MPIerror)
      enddo

C Create new distribution of H and S
      nuog = nuotot
      allocate(Snew(maxnhg))
      call memory('A','D',maxnhg,'diagkp')
      allocate(Hnew(maxnhg,nspin))
      call memory('A','D',maxnhg*nspin,'diagkp')
      allocate(xijloc(3,maxnhg))
      call memory('A','D',3*maxnhg,'diagkp')

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do is = 1,nspin
            do jo = 1,numh(iio)
              Hnew(listhptrg(io)+jo,is) = H(listhptr(iio)+jo,is)
            enddo
          enddo
          do jo = 1,numh(iio)
            Snew(listhptrg(io)+jo) = S(listhptr(iio)+jo)
          enddo
          do jo = 1,numh(iio)
            xijloc(1:3,listhptrg(io)+jo) = xij(1:3,listhptr(iio)+jo)
          enddo
        endif
        do is = 1,nspin
#ifdef NODAT
          call MPI_Bcast(Hnew(listhptrg(io)+1,is),numhg(io),
     .      MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
#else
          call MPI_Bcast(Hnew(listhptrg(io)+1,is),numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)
#endif
        enddo
#ifdef NODAT
        call MPI_Bcast(Snew(listhptrg(io)+1),numhg(io),
     .    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
        call MPI_Bcast(xijloc(1,listhptrg(io)+1),3*numhg(io),
     .    MPI_double_precision,BNode,MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(Snew(listhptrg(io)+1),numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)
        call MPI_Bcast(xijloc(1,listhptrg(io)+1),3*numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)
#endif
      enddo

C Find eigenvalues .........................................
      do ik = 1+Node,nk,Nodes
        do ispin = 1,nspin
          do iuo = 1,nuog
            do juo = 1,nuotot
              Saux(1,juo,iuo) = 0.0d0
              Saux(2,juo,iuo) = 0.0d0
              Haux(1,juo,iuo) = 0.0d0
              Haux(2,juo,iuo) = 0.0d0
            enddo
          enddo
          do io = 1,nuog
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo = listhg(ind)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     .               kpoint(2,ik) * xijloc(2,ind) +
     .               kpoint(3,ik) * xijloc(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
C Note : sign of complex part changed to match change in order of iuo/juo
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+Snew(ind)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-Snew(ind)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+Hnew(ind,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-Hnew(ind,ispin)*skxij
            enddo
          enddo
          call cdiag( Haux, nuotot, Saux, nuotot, nuog,
     .                eo(1,ispin,ik), psi, nuotot, 0, ierror )
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
#ifdef NODAT
        call MPI_Bcast(eo(1,1,ik),nuotot*nspin,MPI_double_precision,
     .    BNode,MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(eo(1,1,ik),nuotot*nspin,DAT_double,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
      enddo

C Check if we are done ................................................
      if (.not.getD) goto 999
C ....................

C Allocate local copies of Dk and Ek for building matrices
      allocate(Dnewloc(maxnhg,2))
      call memory('A','D',2*maxnhg,'diagkp')
      allocate(Enewloc(maxnhg,2))
      call memory('A','D',2*maxnhg,'diagkp')

C Find new Fermi energy and occupation weights ........................
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, maxo, nuotot, eo,
     .               temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, maxspn, nk, wk, maxo, nuotot, eo, 
     .             temp, qtot, qo, ef, Entropy )
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
     .             ( stepf((eo(io,ispin,ik)-e2)/t) -
     .               stepf((eo(io,ispin,ik)-e1)/t)) * 2.0d0 / nspin
            enddo
          enddo
        enddo
      endif
C ....................
      
c New density and energy-density matrices of unit-cell orbitals .......

      Dnewloc(1:maxnhg,1:nspin) = 0.0d0
      Enewloc(1:maxnhg,1:nspin) = 0.0d0

      do ik = 1+Node,nk,Nodes
        do ispin = 1,nspin

C         Find eigenvectors 
          Saux = 0.0d0
          Haux = 0.0d0
          do io = 1,nuog
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo = listhg(ind)
              iuo = indxuo(io)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     .               kpoint(2,ik) * xijloc(2,ind) +
     .               kpoint(3,ik) * xijloc(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Saux(1,juo,iuo)=Saux(1,juo,iuo)+Snew(ind)*ckxij
              Saux(2,juo,iuo)=Saux(2,juo,iuo)-Snew(ind)*skxij
              Haux(1,juo,iuo)=Haux(1,juo,iuo)+Hnew(ind,ispin)*ckxij
              Haux(2,juo,iuo)=Haux(2,juo,iuo)-Hnew(ind,ispin)*skxij
            enddo
          enddo
          call cdiag( Haux, nuotot, Saux, nuotot, nuog,
     .                eo(1,ispin,ik), psi, nuotot, nuog, ierror )
          if (ierror.ne.0) then
            call die('Terminating due to failed diagonalisation')
          endif

C         Add contribution to density matrices of unit-cell orbitals
          Dk = 0.0d0
          Ek = 0.0d0
C Global operation to form new density matrix
          do ie = 1,nuotot
            do j = 1,nuotot
              aux(1,j) = psi(1,j,ie)
              aux(2,j) = psi(2,j,ie)
            enddo
            qe = qo(ie,ispin,ik)
            ee = qo(ie,ispin,ik) * eo(ie,ispin,ik)
            do iuo = 1,nuog
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
          do io = 1,nuog
            iuo = indxuo(io)
C Should be numd/listd/listdptr but H equivalents used to save globalisation
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo = listhg(ind)
              juo = indxuo(jo)
              kxij = kpoint(1,ik) * xijloc(1,ind) +
     .               kpoint(2,ik) * xijloc(2,ind) +
     .               kpoint(3,ik) * xijloc(3,ind)
              ckxij = cos(kxij)
              skxij = sin(kxij)
              Dnewloc(ind,ispin)=Dnewloc(ind,ispin)+ 
     .          Dk(1,juo,iuo)*ckxij - Dk(2,juo,iuo)*skxij
              Enewloc(ind,ispin)=Enewloc(ind,ispin)+ 
     .          Ek(1,juo,iuo)*ckxij - Ek(2,juo,iuo)*skxij
            enddo
          enddo

        enddo
      enddo

C Globalise Dnew/Enew (reference to numhg should be numd)
      do io = 1,nuog
        call WhichNodeOrb(io,Nodes,BNode)
        call GlobalToLocalOrb(io,BNode,Nodes,iio)
        do ispin = 1,nspin
#ifdef NODAT
          call MPI_Reduce(Dnewloc(listhptrg(io)+1,ispin),
     .      Dnew(listhptr(iio)+1,ispin),numhg(io),MPI_double_precision,
     .      MPI_sum,BNode,MPI_Comm_World,MPIerror)
          call MPI_Reduce(Enewloc(listhptrg(io)+1,ispin),
     .      Enew(listhptr(iio)+1,ispin),numhg(io),MPI_double_precision,
     .      MPI_sum,BNode,MPI_Comm_World,MPIerror)
#else
          call MPI_Reduce(Dnewloc(listhptrg(io)+1,ispin),
     .      Dnew(listhptr(iio)+1,ispin),numhg(io),DAT_double,
     .      MPI_sum,BNode,MPI_Comm_World,MPIerror)
          call MPI_Reduce(Enewloc(listhptrg(io)+1,ispin),
     .      Enew(listhptr(iio)+1,ispin),numhg(io),DAT_double,
     .      MPI_sum,BNode,MPI_Comm_World,MPIerror)
#endif
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
      call memory('D','I',size(numhg),'diagkp')
      deallocate(numhg)
      call memory('D','I',size(listhptrg),'diagkp')
      deallocate(listhptrg)
      call memory('D','I',size(listhg),'diagkp')
      deallocate(listhg)

      end
#endif
