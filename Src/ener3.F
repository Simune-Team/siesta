      subroutine ener3(c,grad,lam,eta,enum,h,s,nbasis,nbands,ncmax,
     .                 nctmax,nfmax,maxnh,nhijmax,numc,listc,numct,
     .                 listct,cttoc,numf,listf,numh,listhptr,listh,
     .                 numhij,listhij,ener,nbasisloc)

C ************************************************************************
C Finds the energy at three points of the line passing thru C in the
C direction of GRAD. LAM is the distance (in units of GRAD) between 
C points.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C ****************************** INPUT ***********************************
C real*8 c(ncmax,nbasisloc)    : Current point (wave function coeff.
C                                  in sparse form)
C real*8 grad(ncmax,nbasisloc) : Direction of search (sparse)
C real*8 lam                   : Length of step
C real*8 eta                   : Fermi level parameter of Kim et al.
C real*8 enum                  : Total number of electrons
C real*8 h(maxnh)              : Hamiltonian matrix (sparse)
C real*8 s(maxnh)              : Overlap matrix (sparse)
C integer nbasis               : Number of basis orbitals
C integer nbands               : Number of LWF's
C integer ncmax                : Max num of <>0 elements of each row of C
C integer nctmax               : Max num of <>0 elements of each col of C
C integer nfmax                : Max num of <>0 elements of each row of 
C                                   F = Ct x H
C integer maxnh                : Max num of <>0 elements of H
C integer nhijmax              : Max num of <>0 elements of each row of 
C                                   Hij=Ct x H x C
C integer numc(nbasis)         : Control vector of C matrix
C                                (number of <>0  elements of each row of C)
C integer listc(ncmax,nbasis)  : Control vector of C matrix 
C                               (list of <>0  elements of each row of C)
C integer numct(nbands)        : Control vector of C transpose matrix
C                               (number of <>0  elements of each col of C)
C integer listct(ncmax,nbands) : Control vector of C transpose matrix
C                               (list of <>0  elements of each col of C)
C integer cttoc(ncmax,nbands)  : Map from Ct to C indexing
C integer numf(nbands)         : Control vector of F matrix
C                                (number of <>0  elements of each row of F)
C integer listf(nfmax,nbands)  : Control vector of F matrix 
C                                (list of <>0  elements of each row of F)
C integer numh(nbasisloc)      : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listhptr(nbasisloc)  : Control vector of H matrix
C                                (pointer to the start of each row of H)
C integer listh(maxnh)         : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer numhij(nbands)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C integer nbasisloc            : Number of local atomic orbitals
C ***************************** OUTPUT ***********************************
C real*8 ener(3)               : Energy at the three points:
C                                     C +     lam * GRAD
C                                     C + 2 * lam * GRAD
C                                     C + 3 * lam * GRAD
C ************************************************************************

C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .  nbasis,nbands,ncmax,nctmax,nfmax,maxnh,nhijmax,nbasisloc

      integer
     .  cttoc(nctmax,*),listc(ncmax,nbasis),
     .  listct(nctmax,*),listf(nfmax,*),
     .  listh(maxnh),listhij(nhijmax,*),
     .  numc(nbasis),numct(*),numf(*),
     .  numh(nbasisloc),listhptr(nbasisloc),numhij(*)

      double precision
     .  c(ncmax,nbasisloc),ener(3),eta,enum,grad(ncmax,nbasisloc),
     .  h(maxnh),lam,s(maxnh)

C Internal variables ......................................................

      integer
     .  i,ib,in,ind,j,jn,k,kl,kn,knlh,Nodes,nn,nn2,nbandsmin,nbandsmax,
     .  nbandspernode, iloc, nbandsloc
      integer, dimension(:), allocatable, save :: iptr
      integer, dimension(:), allocatable, save :: iptr2
#ifdef MPI
      integer
     .  MPIerror, nbl, nh, maxnhg, nremainder, ii, Node
      integer, dimension(:), allocatable ::
     .  listhlptr, listhl, numhl
      double precision 
     .  etmp1(6), etmp2(6)
      double precision, dimension(:), allocatable, save :: 
     .  hl,sl
      double precision, dimension(:,:), allocatable, save :: 
     .  cl,gl
      double precision, dimension(:,:), allocatable, save :: 
     .  cl2,gl2
#endif

      double precision
     .  a1,a2,a3,b1,b2,b3,c1,c2,c3,func1(3),func2(3),
     .  lam1,lam2,lam3,pp1,pp2,pp3

      double precision, dimension(:,:), allocatable, save :: 
     .  aux
      double precision, dimension(:,:,:), allocatable, save :: 
     .  bux1, bux2

C..................

#ifdef MPI
C Get Node number
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)

C Get global sizes for first dimension of H/S
      call MPI_AllReduce(maxnh,maxnhg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Allocate local arrays for c/g/h/s
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(cl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'ener3')
      allocate(gl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'ener3')
      allocate(cl2(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'ener3')
      allocate(gl2(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'ener3')
      allocate(hl(maxnhg))
      call memory('A','D',maxnhg,'ener3')
      allocate(sl(maxnhg))
      call memory('A','D',maxnhg,'ener3')

C Initialise hl/sl to prevent possible problems during transfer
      hl(1:maxnhg) = 0.0d0
      sl(1:maxnhg) = 0.0d0

C Work out range of bands for this Node
      nbandspernode = ((nbands-1)/Nodes) 
      nremainder = nbands - nbandspernode*Nodes
      nbandsmin = Node*nbandspernode + min(nremainder,Node) + 1
      nbandsmax = nbandsmin + nbandspernode - 1
      if (Node.lt.nremainder) nbandsmax = nbandsmax + 1
      nbandsloc = nbandsmax - nbandsmin + 1
#else
      Nodes = 1
      nbandsmin = 1
      nbandsmax = nbands
      nbandspernode = nbands
      nbandsloc = nbands
#endif

      call timer('ener3',1)

C Allocate dynamic memory
      allocate(aux(6,nbasis))
      call memory('A','D',6*nbasis,'ener3')
      allocate(bux1(3,nbands,nbandspernode+1))
      call memory('A','D',3*nbands*(nbandspernode+1),'ener3')
      allocate(bux2(3,nbands,nbandspernode+1))
      call memory('A','D',3*nbands*(nbandspernode+1),'ener3')
      allocate(iptr(nbasis))
      call memory('A','I',nbasis,'ener3')
      allocate(iptr2(nbasis))
      call memory('A','I',nbasis,'ener3')

#ifdef MPI
C Allocate local storage for H/S control vectors
      call GetNodeOrbs(nbasis,Node,Nodes,nbl)
      allocate(numhl(nbasis))
      call memory('A','I',nbasis,'ener3')
      allocate(listhlptr(nbasis))
      call memory('A','I',nbasis,'ener3')
      allocate(listhl(maxnhg))
      call memory('A','I',maxnhg,'ener3')
#endif


C Define points to compute energy ..........................................
      lam1 = lam
      lam2 = lam*2.0d0
      lam3 = lam*3.0d0
C..................

C Initialize output and auxiliary varialbles ...............................
      ener(1:3) = 0.0d0
      func1(1:3) = 0.0d0
      func2(1:3) = 0.0d0
      aux(1:6,1:nbasis) = 0.0d0
      bux1(1:3,1:nbands,1:nbandspernode+1) = 0.0d0
      bux2(1:3,1:nbands,1:nbandspernode+1) = 0.0d0

C..................

C Loop over Nodes for broadcasting of c,g,h and s
      do nn = 1,Nodes
        
#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Distribute control vectors for H/S for this node
        if (Node.eq.nn-1) then
          numhl(1:nbl) = numh(1:nbl)
        endif
        call MPI_Bcast(numhl,nbl,MPI_integer,nn-1,
     .    MPI_Comm_World,MPIerror)
        listhlptr(1) = 0
        do in = 2,nbl
          listhlptr(in) = listhlptr(in-1) + numhl(in-1)
        enddo
        nh = listhlptr(nbl) + numhl(nbl)
        if (Node.eq.nn-1) then
          listhl(1:nh) = listh(1:nh)
        endif
        call MPI_Bcast(listhl,nh,MPI_integer,nn-1,
     .    MPI_Comm_World,MPIerror)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          cl(1:ncmax,1:nbl) = c(1:ncmax,1:nbl)
          gl(1:ncmax,1:nbl) = grad(1:ncmax,1:nbl)
          hl(1:nh) = h(1:nh)
          sl(1:nh) = s(1:nh)
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(gl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(hl,nh,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(sl,nh,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

C Build pointer for condensing f and fs on local node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
        enddo

C Loop over Nodes for broadcasting of c,g,h and s
        do nn2 = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
          call GetNodeOrbs(nbasis,nn2-1,Nodes,nbl)

          if (nn.eq.nn2) then
C cl2 and gl2 as needed are already globalised in cl and gl
            cl2(1:ncmax,1:nbl) = cl(1:ncmax,1:nbl)
            gl2(1:ncmax,1:nbl) = gl(1:ncmax,1:nbl)

          else
C Place arrays in local storage on broadcast node
            if (Node.eq.nn2-1) then
              do ib = 1,nbasisloc
                do ii = 1,ncmax
                  cl2(ii,ib) = c(ii,ib)
                  gl2(ii,ib) = grad(ii,ib)
                enddo
              enddo
            endif

C Broadcast arrays if necessary
            call MPI_Bcast(cl2(1,1),ncmax*nbl,DAT_double,nn2-1,
     .          MPI_Comm_World,MPIerror)
            call MPI_Bcast(gl2(1,1),ncmax*nbl,DAT_double,nn2-1,
     .          MPI_Comm_World,MPIerror)
          endif
#endif

C Build pointer for condensing f and fs on local node
          do k = 1,nbasis
            call GlobalToLocalOrb(k,nn2-1,Nodes,kl)
            iptr2(k) = kl
          enddo

C Calculate Functional .....................................................
C F=CtH
C Fs=CtS

          do i = 1,nbandsloc
            do in = 1,numct(i)
              k = listct(in,i)
              kl = iptr(k)
              if (kl.gt.0) then
#ifdef MPI
                pp1 = cl(cttoc(in,i),kl) + lam1 * gl(cttoc(in,i),kl)
                pp2 = cl(cttoc(in,i),kl) + lam2 * gl(cttoc(in,i),kl)
                pp3 = cl(cttoc(in,i),kl) + lam3 * gl(cttoc(in,i),kl)

                do kn = 1,numhl(kl)
                  ind = listhlptr(kl) + kn
                  knlh = listhl(ind)
                  aux(1,knlh) = aux(1,knlh) + pp1 * hl(ind)
                  aux(2,knlh) = aux(2,knlh) + pp2 * hl(ind)
                  aux(3,knlh) = aux(3,knlh) + pp3 * hl(ind)
                  aux(4,knlh) = aux(4,knlh) + pp1 * sl(ind)
                  aux(5,knlh) = aux(5,knlh) + pp2 * sl(ind)
                  aux(6,knlh) = aux(6,knlh) + pp3 * sl(ind)
                enddo
#else
                pp1 = c(cttoc(in,i),kl) + lam1 * grad(cttoc(in,i),kl)
                pp2 = c(cttoc(in,i),kl) + lam2 * grad(cttoc(in,i),kl)
                pp3 = c(cttoc(in,i),kl) + lam3 * grad(cttoc(in,i),kl)
    
                do kn = 1,numh(kl)
                  ind = listhptr(kl) + kn
                  knlh = listh(ind)
                  aux(1,knlh) = aux(1,knlh) + pp1 * h(ind)
                  aux(2,knlh) = aux(2,knlh) + pp2 * h(ind)
                  aux(3,knlh) = aux(3,knlh) + pp3 * h(ind)
                  aux(4,knlh) = aux(4,knlh) + pp1 * s(ind)
                  aux(5,knlh) = aux(5,knlh) + pp2 * s(ind)
                  aux(6,knlh) = aux(6,knlh) + pp3 * s(ind)
                enddo
#endif
              endif
            enddo

            do in = 1,numf(i)
              k = listf(in,i)
              kl = iptr2(k)
              if (kl.gt.0) then
                a1 = aux(1,k)
                a2 = aux(2,k)
                a3 = aux(3,k)
                b1 = aux(4,k)
                b2 = aux(5,k)
                b3 = aux(6,k)

C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC
                do kn = 1,numc(k)
#ifdef MPI
                  c1 = cl2(kn,kl) + lam1 * gl2(kn,kl)
                  c2 = cl2(kn,kl) + lam2 * gl2(kn,kl)
                  c3 = cl2(kn,kl) + lam3 * gl2(kn,kl)
#else
                  c1 = c(kn,kl) + lam1 * grad(kn,kl)
                  c2 = c(kn,kl) + lam2 * grad(kn,kl)
                  c3 = c(kn,kl) + lam3 * grad(kn,kl)
#endif
                  bux1(1,listc(kn,k),i) = bux1(1,listc(kn,k),i) + 
     .              a1 * c1
                  bux1(2,listc(kn,k),i) = bux1(2,listc(kn,k),i) + 
     .              a2 * c2
                  bux1(3,listc(kn,k),i) = bux1(3,listc(kn,k),i) + 
     .              a3 * c3
                  bux2(1,listc(kn,k),i) = bux2(1,listc(kn,k),i) + 
     .              b1 * c1
                  bux2(2,listc(kn,k),i) = bux2(2,listc(kn,k),i) + 
     .              b2 * c2
                  bux2(3,listc(kn,k),i) = bux2(3,listc(kn,k),i) + 
     .              b3 * c3
                enddo
              endif
            enddo

            do in = 1,numct(i)
              k = listct(in,i)
              kl = iptr(k)
              if (kl.gt.0) then
#ifdef MPI
                do kn = 1,numhl(kl)
                  ind = listhlptr(kl) + kn
                  knlh = listhl(ind)
                  aux(1,knlh) = 0.0d0
                  aux(2,knlh) = 0.0d0
                  aux(3,knlh) = 0.0d0
                  aux(4,knlh) = 0.0d0
                  aux(5,knlh) = 0.0d0
                  aux(6,knlh) = 0.0d0
                enddo
#else
                do kn = 1,numh(kl)
                  ind = listhptr(kl) + kn
                  knlh = listh(ind)
                  aux(1,knlh) = 0.0d0
                  aux(2,knlh) = 0.0d0
                  aux(3,knlh) = 0.0d0
                  aux(4,knlh) = 0.0d0
                  aux(5,knlh) = 0.0d0
                  aux(6,knlh) = 0.0d0
                enddo
#endif
              endif
            enddo

C End of loop over nbandsloc
          enddo

C End of loop over node nn2
        enddo

C End of loop over node nn
      enddo

      do i = 1,nbandsloc
        iloc = i + nbandsmin - 1
        func1(1) = func1(1) + bux1(1,iloc,i) - eta * bux2(1,iloc,i)
        func1(2) = func1(2) + bux1(2,iloc,i) - eta * bux2(2,iloc,i)
        func1(3) = func1(3) + bux1(3,iloc,i) - eta * bux2(3,iloc,i)
        do jn = 1,numhij(i)
          func2(1) = func2(1) + (bux1(1,listhij(jn,i),i) -
     .             eta * bux2(1,listhij(jn,i),i)) * 
     .             bux2(1,listhij(jn,i),i)
          func2(2) = func2(2) + (bux1(2,listhij(jn,i),i) -
     .             eta * bux2(2,listhij(jn,i),i)) * 
     .             bux2(2,listhij(jn,i),i)
          func2(3) = func2(3) + (bux1(3,listhij(jn,i),i) -
     .             eta * bux2(3,listhij(jn,i),i)) * 
     .             bux2(3,listhij(jn,i),i)
        enddo
      enddo


#ifdef MPI
C Globalisation of func1 and func2
      do i=1,3
        etmp1(i) = func1(i)
        etmp1(3+i) = func2(i)
      enddo
      call MPI_AllReduce(etmp1,etmp2,6,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      do i=1,3
        func1(i) = etmp2(i)
        func2(i) = etmp2(3+i)
      enddo
#endif

C This is valid for an spin-unpolarized sytem
      do i=1,3
        ener(i) = 2.0d0 * func1(i) - func2(i) + eta * enum / 2.0d0
      enddo

C ...................

C Free dynamic memory
#ifdef MPI
      call memory('D','D',size(cl),'ener3')
      deallocate(cl)
      call memory('D','D',size(gl),'ener3')
      deallocate(gl)
      call memory('D','D',size(cl2),'ener3')
      deallocate(cl2)
      call memory('D','D',size(gl2),'ener3')
      deallocate(gl2)
      call memory('D','D',size(hl),'ener3')
      deallocate(hl)
      call memory('D','D',size(sl),'ener3')
      deallocate(sl)
      call memory('D','I',size(numhl),'ener3')
      deallocate(numhl)
      call memory('D','I',size(listhlptr),'ener3')
      deallocate(listhlptr)
      call memory('D','I',size(listhl),'ener3')
      deallocate(listhl)
#endif
      call memory('D','D',size(aux),'ener3')
      deallocate(aux)
      call memory('D','D',size(bux1),'ener3')
      deallocate(bux1)
      call memory('D','D',size(bux2),'ener3')
      deallocate(bux2)
      call memory('D','I',size(iptr),'ener3')
      deallocate(iptr)
      call memory('D','I',size(iptr2),'ener3')
      deallocate(iptr2)

      call timer('ener3',2)
      return
      end
