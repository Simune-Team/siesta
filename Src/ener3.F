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
     .  i,ib,in,ind,j,jn,k,kl,kn,knk,knlh,Nodes,nn,nbandsmin,
     .  nbandsmax,nbandspernode,iloc,nbandsloc,nfmaxg,nfmaxl,
     .  nhijmaxg
      integer, dimension(:), allocatable, save :: 
     .  iptr
#ifdef MPI
      integer
     .  MPIerror, nbl, nh, maxnhg, nhijmaxl, nremainder, ii, Node
      integer, dimension(:), allocatable ::
     .  listhlptr, listhl, numhl
      double precision 
     .  etmp1(6), etmp2(6)
      double precision, dimension(:), allocatable, save :: 
     .  hl,sl
      double precision, dimension(:,:), allocatable, save :: 
     .  cl,gl
#endif

      double precision
     .  a1,a2,a3,b1,b2,b3,c1,c2,c3,func1(3),func2(3),
     .  lam1,lam2,lam3,pp1,pp2,pp3

      double precision, dimension(:,:), allocatable, save :: 
     .  dux, cux1, cux2
      double precision, dimension(:,:,:), allocatable, save :: 
     .  aux, bux1, bux2

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

C Find true value of nfmax
      nfmaxl = 0
      do i = 1,nbandsloc
        nfmaxl = max(nfmaxl,numf(i))
      enddo
      call MPI_AllReduce(nfmaxl,nfmaxg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Find true value of nhijmax
      nhijmaxl = 0
      do i = 1,nbandsloc
        nhijmaxl = max(nhijmaxl,numhij(i))
      enddo
      call MPI_AllReduce(nhijmaxl,nhijmaxg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)
#else
      Nodes = 1
      nbandsmin = 1
      nbandsmax = nbands
      nbandspernode = nbands
      nbandsloc = nbands
      nfmaxg = nfmax
      nhijmaxg = nhijmax
#endif

      call timer('ener3',1)

C Allocate dynamic memory
      allocate(aux(6,nfmaxg,nbandsloc))
      call memory('A','D',6*nfmaxg*nbandsloc,'ener3')
      allocate(dux(6,nbasis))
      call memory('A','D',6*nbasis,'ener3')
      allocate(iptr(nbasis))
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
      dux(1:6,1:nbasis) = 0.0d0
      aux(1:6,1:nfmaxg,1:nbandsloc) = 0.0d0

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

C Calculate Functional .....................................................
C F=CtH
C Fs=CtS

        do i = 1,nbandsloc
          iloc = i + nbandsmin - 1
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
                dux(1,knlh) = dux(1,knlh) + pp1 * hl(ind)
                dux(2,knlh) = dux(2,knlh) + pp2 * hl(ind)
                dux(3,knlh) = dux(3,knlh) + pp3 * hl(ind)
                dux(4,knlh) = dux(4,knlh) + pp1 * sl(ind)
                dux(5,knlh) = dux(5,knlh) + pp2 * sl(ind)
                dux(6,knlh) = dux(6,knlh) + pp3 * sl(ind)
              enddo
#else
              pp1 = c(cttoc(in,i),kl) + lam1 * grad(cttoc(in,i),kl)
              pp2 = c(cttoc(in,i),kl) + lam2 * grad(cttoc(in,i),kl)
              pp3 = c(cttoc(in,i),kl) + lam3 * grad(cttoc(in,i),kl)
    
              do kn = 1,numh(kl)
                ind = listhptr(kl) + kn
                knlh = listh(ind)
                dux(1,knlh) = dux(1,knlh) + pp1 * h(ind)
                dux(2,knlh) = dux(2,knlh) + pp2 * h(ind)
                dux(3,knlh) = dux(3,knlh) + pp3 * h(ind)
                dux(4,knlh) = dux(4,knlh) + pp1 * s(ind)
                dux(5,knlh) = dux(5,knlh) + pp2 * s(ind)
                dux(6,knlh) = dux(6,knlh) + pp3 * s(ind)
              enddo
#endif
            endif
          enddo

C Reduce dux to elements needed in aux
          do in = 1,numf(i)
            k = listf(in,i)
            aux(1,in,i) = aux(1,in,i) + dux(1,k)
            aux(2,in,i) = aux(2,in,i) + dux(2,k)
            aux(3,in,i) = aux(3,in,i) + dux(3,k)
            aux(4,in,i) = aux(4,in,i) + dux(4,k)
            aux(5,in,i) = aux(5,in,i) + dux(5,k)
            aux(6,in,i) = aux(6,in,i) + dux(6,k)
          enddo

          do in = 1,numct(i)
            k = listct(in,i)
            kl = iptr(k)
            if (kl.gt.0) then
#ifdef MPI
              do kn = 1,numhl(kl)
                ind = listhlptr(kl) + kn
                knlh = listhl(ind)
                dux(1,knlh) = 0.0d0
                dux(2,knlh) = 0.0d0
                dux(3,knlh) = 0.0d0
                dux(4,knlh) = 0.0d0
                dux(5,knlh) = 0.0d0
                dux(6,knlh) = 0.0d0
              enddo
#else
              do kn = 1,numh(kl)
                ind = listhptr(kl) + kn
                knlh = listh(ind)
                dux(1,knlh) = 0.0d0
                dux(2,knlh) = 0.0d0
                dux(3,knlh) = 0.0d0
                dux(4,knlh) = 0.0d0
                dux(5,knlh) = 0.0d0
                dux(6,knlh) = 0.0d0
              enddo
#endif
            endif
          enddo

C End of loop over nbandsloc
        enddo

C End of loop over node nn
      enddo

C Free local arrays that are not needed further
      call memory('D','D',size(dux),'ener3')
      deallocate(dux)
#ifdef MPI
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

C Allocate and initialise local arrays cux1/cux2 
      allocate(cux1(3,nbands))
      call memory('A','D',3*nbands,'ener3')
      allocate(cux2(3,nbands))
      call memory('A','D',3*nbands,'ener3')
      cux1(1:3,1:nbands) = 0.0d0
      cux2(1:3,1:nbands) = 0.0d0

C Allocate and initialise bux1/bux2
      allocate(bux1(3,nhijmaxg,nbandspernode+1))
      call memory('A','D',3*nhijmaxg*(nbandspernode+1),'ener3')
      allocate(bux2(3,nhijmaxg,nbandspernode+1))
      call memory('A','D',3*nhijmaxg*(nbandspernode+1),'ener3')
      bux1(1:3,1:nhijmaxg,1:nbandspernode+1) = 0.0d0
      bux2(1:3,1:nhijmaxg,1:nbandspernode+1) = 0.0d0

C Loop over Nodes for broadcasting of c and g
      do nn = 1,Nodes
        
#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do ib = 1,nbl
            do ii = 1,ncmax
              cl(ii,ib) = c(ii,ib)
              gl(ii,ib) = grad(ii,ib)
            enddo
          enddo
        endif

C Broadcast arrays if necessary
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .      MPI_Comm_World,MPIerror)
        call MPI_Bcast(gl(1,1),ncmax*nbl,DAT_double,nn-1,
     .      MPI_Comm_World,MPIerror)
#endif

C Build pointer for condensing f and fs on local node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
        enddo

C Calculate Functional .....................................................
C F=CtH
C Fs=CtS

        do i = 1,nbandsloc
          iloc = i + nbandsmin - 1

          do in = 1,numf(i)
            k = listf(in,i)
            kl = iptr(k)
            if (kl.gt.0) then
              a1 = aux(1,in,i)
              a2 = aux(2,in,i)
              a3 = aux(3,in,i)
              b1 = aux(4,in,i)
              b2 = aux(5,in,i)
              b3 = aux(6,in,i)

C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC
              do kn = 1,numc(k)
                knk = listc(kn,k)
#ifdef MPI
                c1 = cl(kn,kl) + lam1 * gl(kn,kl)
                c2 = cl(kn,kl) + lam2 * gl(kn,kl)
                c3 = cl(kn,kl) + lam3 * gl(kn,kl)
#else
                c1 = c(kn,kl) + lam1 * grad(kn,kl)
                c2 = c(kn,kl) + lam2 * grad(kn,kl)
                c3 = c(kn,kl) + lam3 * grad(kn,kl)
#endif
                cux1(1,knk) = cux1(1,knk) + a1 * c1
                cux1(2,knk) = cux1(2,knk) + a2 * c2
                cux1(3,knk) = cux1(3,knk) + a3 * c3
                cux2(1,knk) = cux2(1,knk) + b1 * c1
                cux2(2,knk) = cux2(2,knk) + b2 * c2
                cux2(3,knk) = cux2(3,knk) + b3 * c3
              enddo
            endif
          enddo

C Add diagonal term to energy
          func1(1) = func1(1) + cux1(1,iloc) - eta * cux2(1,iloc)
          func1(2) = func1(2) + cux1(2,iloc) - eta * cux2(2,iloc)
          func1(3) = func1(3) + cux1(3,iloc) - eta * cux2(3,iloc)

C Add cux to bux arrays and re-zero cux
          do jn = 1,numhij(i)
            j = listhij(jn,i)
            bux1(1,jn,i) = bux1(1,jn,i) + cux1(1,j)
            bux1(2,jn,i) = bux1(2,jn,i) + cux1(2,j)
            bux1(3,jn,i) = bux1(3,jn,i) + cux1(3,j)
            bux2(1,jn,i) = bux2(1,jn,i) + cux2(1,j)
            bux2(2,jn,i) = bux2(2,jn,i) + cux2(2,j)
            bux2(3,jn,i) = bux2(3,jn,i) + cux2(3,j)
            cux1(1,j) = 0.0d0
            cux1(2,j) = 0.0d0
            cux1(3,j) = 0.0d0
            cux2(1,j) = 0.0d0
            cux2(2,j) = 0.0d0
            cux2(3,j) = 0.0d0
          enddo

C End of loop over nbandsloc
        enddo

C End of loop over node nn
      enddo

      do i = 1,nbandsloc
        do jn = 1,numhij(i)
          func2(1) = func2(1) + (bux1(1,jn,i) - eta * bux2(1,jn,i)) * 
     .             bux2(1,jn,i)
          func2(2) = func2(2) + (bux1(2,jn,i) - eta * bux2(2,jn,i)) * 
     .             bux2(2,jn,i)
          func2(3) = func2(3) + (bux1(3,jn,i) - eta * bux2(3,jn,i)) * 
     .             bux2(3,jn,i)
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
#endif
      call memory('D','D',size(aux),'ener3')
      deallocate(aux)
      call memory('D','D',size(cux1),'ener3')
      deallocate(cux1)
      call memory('D','D',size(cux2),'ener3')
      deallocate(cux2)
      call memory('D','D',size(bux1),'ener3')
      deallocate(bux1)
      call memory('D','D',size(bux2),'ener3')
      deallocate(bux2)
      call memory('D','I',size(iptr),'ener3')
      deallocate(iptr)

      call timer('ener3',2)
      return
      end
