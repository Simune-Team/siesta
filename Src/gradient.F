      subroutine gradient(c,eta,enum,h,s,nbasis,nbands,ncmax,nctmax,
     .                    nfmax,nftmax,maxnh,nhijmax,numc,listc,
     .                    numct,listct,cttoc,numf,listf,numft,listft,
     .                    fttof,numh,listhptr,listh,numhij,listhij,
     .                    indgf,grad,ener,nbasisloc)

C ************************************************************************
C Finds the energy and gradient at point C.
C Uses the functional of Kim et al (PRB 52, 1640 (95))
C Works only with spin-unpolarized systems
C Written by P.Ordejon. October'96
C Last modified: J.M.Soler. 30/04/97
C Parallelisation introduced by J.D. Gale, April'99
C ****************************** INPUT ***********************************
C real*8 c(ncmax,nbasisloc)    : Current point (wave function coeff.
C                                  in sparse form)
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
C integer nftmax               : Max num of <>0 elements of each col of F
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
C integer numft(nbasis)        : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,nbasisloc): Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer fttof(nfmax,nbasisloc): Map from Ft to F indexing
C integer numh(nbasisloc)      : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listhptr(nbasisloc)  : Control vector of H matrix
C                                (pointer to the start of each row of H)
C integer listh(maxnh)         : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer numhij(nbasis)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbandsloc):Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
C integer indgf(ncmax,nbasisloc): Maps elements of GRAD and F
C ***************************** OUTPUT ***********************************
C real*8 ener                  : Energy at point C
C real*8 grad(ncmax,nbasisloc) : Gradient of functional (sparse)
C ************************************************************************
C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer
     .  nbasis,nbands,ncmax,nctmax,nfmax,nftmax,maxnh,nhijmax,nbasisloc

      integer
     .  cttoc(nctmax,*),indgf(ncmax,nbasisloc),
     .  listc(ncmax,nbasis),listct(nctmax,*),listf(nfmax,*),
     .  listft(nftmax,nbasisloc),listh(maxnh),listhij(nhijmax,*),
     .  numc(nbasis),numct(*),numf(*),numft(nbasis),numh(nbasisloc),
     .  listhptr(nbasisloc),numhij(*),fttof(nftmax,nbasisloc)

      double precision
     .  c(ncmax,nbasisloc),ener,eta,enum,grad(ncmax,nbasisloc),
     .  h(maxnh),s(maxnh)

C Internal variables ......................................................

      integer
     .  i,ik,in,ind,j,jk,jn,k,kl,kn,knk,knlh,mu,muk,Node,Nodes,
     .  nn, nbandsmin, nbandsmax, nbandspernode, iloc, 
     .  nbandsloc, nfmaxg, nhijmaxg
      integer, dimension(:), allocatable, save :: 
     .  iptr, iptrlg
#ifdef MPI
      integer
     .  MPIerror, ii, il, jloc, mul, nh, nremainder, 
     .  nbmin, nbmax, maxnhg, nloc, nbl, nfmaxl,
     .  nhijmaxl
      integer, dimension(:), allocatable ::
     .  listhlptr, listhl, numhl
      integer, dimension(:,:), allocatable, save :: 
     .  listftl, indgfl
      double precision
     .  rtmp1(2),rtmp2(2)
      double precision, dimension(:), allocatable, save :: 
     .  hl, sl
      double precision, dimension(:,:), allocatable, save :: 
     .  cl, gl, ftl, ftl2
#endif
      double precision, dimension(:), allocatable, save :: 
     .  aux1, aux2
      double precision, dimension(:,:), allocatable, save :: 
     .  f, ft, fs, fst, bux1, bux2

      double precision
     .  a0,b0,p0,func1,func2
     
C..................

      call timer('gradient',1)
      
C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)

C Get global sizes for first dimension of H/S
      call MPI_AllReduce(maxnh,maxnhg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Allocate local arrays for c/h/s
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(cl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'gradient')
      allocate(hl(maxnhg))
      call memory('A','D',maxnhg,'gradient')
      allocate(sl(maxnhg))
      call memory('A','D',maxnhg,'gradient')

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

C Get global maximum number of non-zero elements in a given row
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
      Node = 0
      Nodes = 1
      nbandsmin = 1
      nbandsmax = nbands
      nbandspernode = nbands
      nbandsloc = nbands
      nfmaxg = nfmax
      nhijmaxg = nhijmax
#endif

C Allocate aux arrays
      allocate(aux1(nbasis))
      call memory('A','D',nbasis,'gradient')
      allocate(aux2(nbasis))
      call memory('A','D',nbasis,'gradient')
      aux1(1:nbasis) = 0.0d0
      aux2(1:nbasis) = 0.0d0

C Allocate dynamic memory
      allocate(iptr(nbasis))
      call memory('A','I',nbasis,'gradient')
      allocate(iptrlg(nbasis))
      call memory('A','I',nbasisloc,'gradient')

C Setup local to global pointer
      do i = 1,nbasisloc
        call LocalToGlobalOrb(i,Node,Nodes,iptrlg(i))
      enddo

#ifdef MPI
C Allocate local storage for H/S control vectors
      call GetNodeOrbs(nbasis,Node,Nodes,nbl)
      allocate(numhl(nbasis))
      call memory('A','I',nbasis,'gradient')
      allocate(listhlptr(nbasis))
      call memory('A','I',nbasis,'gradient')
      allocate(listhl(maxnhg))
      call memory('A','I',maxnhg,'gradient')
#endif

C Initialize output and auxiliary varialbles ...............................

      ener = 0.0d0
      func1 = 0.0d0
      func2 = 0.0d0

      allocate(f(nfmaxg,nbandspernode+1))
      call memory('A','D',nfmaxg*(nbandspernode+1),'gradient')
      allocate(fs(nfmaxg,nbandspernode+1))
      call memory('A','D',nfmaxg*(nbandspernode+1),'gradient')
      allocate(bux1(nhijmaxg,nbandspernode+1))
      call memory('A','D',nhijmaxg*(nbandspernode+1),'gradient')
      allocate(bux2(nhijmaxg,nbandspernode+1))
      call memory('A','D',nhijmaxg*(nbandspernode+1),'gradient')

      f(1:nfmaxg,1:nbandspernode+1) = 0.0d0
      fs(1:nfmaxg,1:nbandspernode+1) = 0.0d0
      bux1(1:nhijmaxg,1:nbandspernode+1) = 0.0d0
      bux2(1:nhijmaxg,1:nbandspernode+1) = 0.0d0

C..................

C Calculate Functional .....................................................

C Start loop over Nodes
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
          hl(1:nh) = h(1:nh)
          sl(1:nh) = s(1:nh)
        endif

C Broadcast arrays
#ifdef NODAT
        call MPI_Bcast(cl(1,1),ncmax*nbl,MPI_double_precision,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(hl,maxnhg,MPI_double_precision,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(sl,maxnhg,MPI_double_precision,nn-1,
     .        MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(hl,maxnhg,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(sl,maxnhg,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

C Build pointer to local orbitals for current working Node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,iptr(k))
        enddo
#endif

C F=CtH  ---> JMS: F=Ct*(H-eta*S)
C Fs=CtS
        do i = 1,nbandsloc
          do in = 1,numct(i)
            k = listct(in,i)
#ifdef MPI
            kl = iptr(k)
            if (kl.gt.0) then
              ik = cttoc(in,i)
              p0 = cl(ik,kl)
              ind = listhlptr(kl)
              do kn = 1,numhl(kl)
                ind = ind + 1
                aux1(listhl(ind)) = aux1(listhl(ind)) +
     .            p0 * (hl(ind)-eta*sl(ind))
                aux2(listhl(ind)) = aux2(listhl(ind)) +
     .            p0 * sl(ind)
              enddo
            endif
#else
            ik = cttoc(in,i)
            p0 = c(ik,k)
            ind = listhptr(k)
            do kn = 1,numh(k)
              ind = ind + 1
              aux1(listh(ind)) = aux1(listh(ind)) +
     .          p0 * (h(ind)-eta*s(ind))
              aux2(listh(ind)) = aux2(listh(ind)) +
     .          p0 * s(ind)
            enddo
#endif
          enddo
          do in = 1,numf(i)
            k = listf(in,i)
            f(in,i) = f(in,i) + aux1(k)
            fs(in,i) = fs(in,i) + aux2(k)
            aux1(k) = 0.0d0
            aux2(k) = 0.0d0
          enddo
        enddo

C End loop over Nodes
      enddo

#ifdef MPI
C Deallocate memory for hl and sl as no longer needed
      call memory('D','D',size(hl),'gradient')
      deallocate(hl)
      call memory('D','D',size(sl),'gradient')
      deallocate(sl)
      call memory('D','I',size(numhl),'gradient')
      deallocate(numhl)
      call memory('D','I',size(listhlptr),'gradient')
      deallocate(listhlptr)
      call memory('D','I',size(listhl),'gradient')
      deallocate(listhl)
#endif

C Start loop over Nodes
      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          cl(1:ncmax,1:nbl) = c(1:ncmax,1:nbl)
        endif

C Broadcast arrays
#ifdef NODAT
        call MPI_Bcast(cl(1,1),ncmax*nbl,MPI_double_precision,nn-1,
     .        MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

C Build pointer to local orbitals for current working Node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
        enddo
#endif

C Hij=CtHC
C Sij=CtSC
C multiply FxC and FsxC row by row

        do i = 1,nbandsloc

          iloc = i + nbandsmin - 1

          do in = 1,numf(i)
            k = listf(in,i)
#ifdef MPI
            kl = iptr(k)
            if (kl.gt.0) then
#endif
              a0 = f(in,i)
              b0 = fs(in,i)
              do kn = 1,numc(k)
                knk = listc(kn,k)
#ifdef MPI
                aux1(knk) = aux1(knk) + a0 * cl(kn,kl)
                aux2(knk) = aux2(knk) + b0 * cl(kn,kl)
#else
                aux1(knk) = aux1(knk) + a0 * c(kn,k)
                aux2(knk) = aux2(knk) + b0 * c(kn,k)
#endif
              enddo
#ifdef MPI
            endif
#endif
          enddo

C Evaluate first energy terms
          func1 = func1 + aux1(iloc)

C Condense aux down to bux and zero aux again
          do jn = 1,numhij(i)
            j = listhij(jn,i)
            bux1(jn,i) = bux1(jn,i) + aux1(j)
            bux2(jn,i) = bux2(jn,i) + aux2(j)
            aux1(j) = 0.0d0
            aux2(j) = 0.0d0
          enddo

        enddo

C End loop over Nodes
      enddo

C Calculate second part of energy
      do i = 1,nbandsloc
        do jn = 1,numhij(i)
          func2 = func2 + bux1(jn,i) * bux2(jn,i)
        enddo
      enddo

#ifdef MPI
C Allocate local storage for broadcast copy of indgf
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(indgfl(ncmax,nbl))
      call memory('A','I',ncmax*nbl,'gradient')

C Add contribution from f to grad
      do nn = 1,Nodes

C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)
        cl(1:ncmax,1:nbl) = 0.0

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          indgfl(1:ncmax,1:nbl) = indgf(1:ncmax,1:nbl)
        endif

C Broadcast arrays
        call MPI_Bcast(indgfl(1,1),ncmax*nbl,MPI_integer,nn-1,
     .    MPI_Comm_World,MPIerror)

        do kl = 1,nbl
          call LocalToGlobalOrb(kl,nn-1,Nodes,k)
          do ik = 1,numc(k)
            i = indgfl(ik,kl)
            j = listc(ik,k)
            if (j.ge.nbandsmin.and.j.le.nbandsmax) then
              jloc = j - nbandsmin + 1
              cl(ik,kl) = 4.d0 * f(i,jloc) + cl(ik,kl)
            endif
          enddo
        enddo

C Global reduction of grad
#ifdef NODAT
        call MPI_Reduce(cl(1,1),grad(1,1),ncmax*nbl,
     .    MPI_double_precision,MPI_sum,nn-1,
     .    MPI_Comm_World,MPIerror)
#else
        call MPI_Reduce(cl(1,1),grad(1,1),ncmax*nbl,
     .    DAT_double,MPI_sum,nn-1,
     .    MPI_Comm_World,MPIerror)
#endif

      enddo

C Free local storage for indgf
      call memory('D','I',size(indgfl),'gradient')
      deallocate(indgfl)
#else
      do k = 1,nbasis
        do ik = 1,numc(k)
          i = indgf(ik,k)
          j = listc(ik,k)
          grad(ik,k) = 4.d0 * f(i,j)
        enddo
      enddo
#endif

C Build pointer to local orbitals for current working Node
      do k = 1,nbasis
        call GlobalToLocalOrb(k,Node,Nodes,kl)
        iptr(k) = kl
      enddo

C Allocate memory for transpose arrays ft and fst
      allocate(ft(nftmax,nbasisloc))
      call memory('A','D',nftmax*nbasisloc,'gradient')
#ifdef MPI
      allocate(ftl(nfmaxg,nbandspernode+1))
      call memory('A','D',nfmaxg*(nbandspernode+1),'gradient')
#endif

C Calculate transpose of f and fs
#ifdef MPI
      do nn = 1,Nodes

C Get local size on broadcast node
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nloc = (nbmax-nbmin+1)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do i = 1,nloc
            do kn = 1,numf(i)
              ftl(kn,i) = f(kn,i)
            enddo
          enddo
        endif

C Broadcast arrays
#ifdef NODAT
        call MPI_Bcast(ftl(1,1),nfmaxg*nloc,MPI_double_precision,nn-1,
     .    MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(ftl(1,1),nfmaxg*nloc,DAT_double,nn-1,
     .    MPI_Comm_World,MPIerror)
#endif

        do mu = 1,nbasisloc
          do muk = 1,numft(iptrlg(mu))
            j = listft(muk,mu)
            if (j.ge.nbmin.and.j.le.nbmax) then
              iloc = j - nbmin + 1
              jk = fttof(muk,mu)
              ft(muk,mu) = ftl(jk,iloc)
            endif
          enddo
        enddo
      enddo
#else
      do mu = 1,nbasis
        do muk = 1,numft(mu)
          j = listft(muk,mu)
          jk = fttof(muk,mu)
          ft(muk,mu) = f(jk,j)
        enddo
      enddo
#endif

C Deallocate f
      call memory('D','D',size(f),'gradient')
      deallocate(f)

C Allocate transpose of fs
      allocate(fst(nftmax,nbasisloc))
      call memory('A','D',nftmax*nbasisloc,'gradient')

#ifdef MPI
      do nn = 1,Nodes

C Get local size on broadcast node
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nloc = (nbmax-nbmin+1)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do i = 1,nloc
            do kn = 1,numf(i)
              ftl(kn,i) = fs(kn,i)
            enddo
          enddo
        endif

C Broadcast arrays
#ifdef NODAT
        call MPI_Bcast(ftl(1,1),nfmaxg*nloc,MPI_double_precision,nn-1,
     .    MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(ftl(1,1),nfmaxg*nloc,DAT_double,nn-1,
     .    MPI_Comm_World,MPIerror)
#endif

        do mu = 1,nbasisloc
          do muk = 1,numft(iptrlg(mu))
            j = listft(muk,mu)
            if (j.ge.nbmin.and.j.le.nbmax) then
              iloc = j - nbmin + 1
              jk = fttof(muk,mu)
              fst(muk,mu) = ftl(jk,iloc)
            endif
          enddo
        enddo
      enddo

C Free ftl with old dimensions 
      call memory('D','D',size(ftl),'gradient')
      deallocate(ftl)
#else
      do mu = 1,nbasis
        do muk = 1,numft(mu)
          j = listft(muk,mu)
          jk = fttof(muk,mu)
          fst(muk,mu) = fs(jk,j)
        enddo
      enddo
#endif

C Deallocate fs as it is not used again
      call memory('D','D',size(fs),'gradient')
      deallocate(fs)

#ifdef MPI
C Create local storage memory for transpose arrays
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(ftl2(nftmax,nbl))
      call memory('A','D',nftmax*nbl,'gradient')
      allocate(listftl(nftmax,nbl))
      call memory('A','I',nftmax*nbl,'gradient')
      allocate(gl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'gradient')
#endif

C Re-zero aux arrays
      aux1(1:nbasis) = 0.0d0
      aux2(1:nbasis) = 0.0d0

      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)
        cl(1:ncmax,1:nbl) = 0.0

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          ftl2(1:nftmax,1:nbasisloc) = ft(1:nftmax,1:nbasisloc)
          do i = 1,nbasisloc
            call LocalToGlobalOrb(i,Node,Nodes,il)
            do ii = 1,numft(il)
              listftl(ii,i) = listft(ii,i)
            enddo
            do ii = numft(il)+1,nftmax
              listftl(ii,i) = 0
            enddo
          enddo
        endif

C Broadcast arrays
#ifdef NODAT
        call MPI_Bcast(ftl2(1,1),nftmax*nbl,MPI_double_precision,nn-1,
     .        MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(ftl2(1,1),nftmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif
        call MPI_Bcast(listftl(1,1),nftmax*nbl,MPI_integer,
     .        nn-1,MPI_Comm_World,MPIerror)

C Build pointer to local orbitals for current working Node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
        enddo
#endif
        do i = 1,nbandsloc

C Expand bux into aux
          do jn = 1,numhij(i)
            j = listhij(jn,i)
            aux1(j) = bux1(jn,i)
            aux2(j) = bux2(jn,i)
          enddo

C multiply Hij x Fs and Sij x F row by row
C (only products of neccesary elements)
#ifdef MPI
          do ik = 1,numct(i)
            mu = listct(ik,i)
            mul = iptr(mu)
            if (mul.gt.0) then
              a0 = 0.d0
              do muk = 1,numft(mu)
                j = listftl(muk,mul)
                a0 = a0 + aux2(j) * ftl2(muk,mul)
              enddo
              cl(cttoc(ik,i),mul) = cl(cttoc(ik,i),mul) -2.0d0*a0
            endif
          enddo
#else
          do ik = 1,numct(i)
            mu = listct(ik,i)
            a0 = 0.d0
            do muk = 1,numft(mu)
              j = listft(muk,mu)
              a0 = a0 + aux1(j) * fst(muk,mu) + aux2(j) * ft(muk,mu)
            enddo
            grad(cttoc(ik,i),mu) = grad(cttoc(ik,i),mu) -2.0d0*a0
          enddo
#endif

C Re-zero aux arrays
          do jn = 1,numhij(i)
            j = listhij(jn,i)
            aux1(j) = 0.0d0
            aux2(j) = 0.0d0
          enddo

        enddo

#ifdef MPI
C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          ftl2(1:nftmax,1:nbl) = fst(1:nftmax,1:nbl)
        endif

C Broadcast arrays
#ifdef NODAT
        call MPI_Bcast(ftl2(1,1),nftmax*nbl,MPI_double_precision,nn-1,
     .        MPI_Comm_World,MPIerror)
#else
        call MPI_Bcast(ftl2(1,1),nftmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

        do i = 1,nbandsloc

C Expand bux into aux
          do jn = 1,numhij(i)
            j = listhij(jn,i)
            aux1(j) = bux1(jn,i)
          enddo

C multiply Hij x Fs and Sij x F row by row
C (only products of neccesary elements)
          do ik = 1,numct(i)
            mu = listct(ik,i)
            mul = iptr(mu)
            if (mul.gt.0) then
              a0 = 0.d0
              do muk = 1,numft(mu)
                j = listftl(muk,mul)
                a0 = a0 + aux1(j) * ftl2(muk,mul)
              enddo
              cl(cttoc(ik,i),mul) = cl(cttoc(ik,i),mul) -2.0 * a0
            endif
          enddo

C Rezero aux1
          do jn = 1,numhij(i)
            j = listhij(jn,i)
            aux1(j) = 0.0d0
          enddo

        enddo
#endif

#ifdef MPI
C Global reduction of grad
#ifdef NODAT
        call MPI_Reduce(cl(1,1),gl(1,1),ncmax*nbl,
     .    MPI_double_precision,MPI_sum,nn-1,
     .    MPI_Comm_World,MPIerror)
#else
        call MPI_Reduce(cl(1,1),gl(1,1),ncmax*nbl,
     .    DAT_double,MPI_sum,nn-1,
     .    MPI_Comm_World,MPIerror)
#endif
        if (Node.eq.nn-1) then
          do kl = 1,nbasisloc
            call LocalToGlobalOrb(kl,nn-1,Nodes,k)
            do ik = 1,numc(k)
              grad(ik,kl) = gl(ik,kl) + grad(ik,kl)
            enddo
          enddo
        endif
#endif

C End loop over Nodes
      enddo

#ifdef MPI
C Global reduction of func1 and func2
      rtmp1(1)=func1
      rtmp1(2)=func2
#ifdef NODAT
      call MPI_AllReduce(rtmp1,rtmp2,2,MPI_double_precision,MPI_sum,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(rtmp1,rtmp2,2,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
#endif
      func1 = rtmp2(1)
      func2 = rtmp2(2)
#endif
      ener = 2.d0 * func1 - func2 + eta * enum / 2.d0
C ...................

C Deallocate aux1/aux2
      call memory('D','D',size(aux1),'gradient')
      deallocate(aux1)
      call memory('D','D',size(aux2),'gradient')
      deallocate(aux2)

C Free dynamic memory
      call memory('D','I',size(iptr),'gradient')
      deallocate(iptr)
      call memory('D','I',size(iptrlg),'gradient')
      deallocate(iptrlg)
      call memory('D','D',size(ft),'gradient')
      deallocate(ft)
      call memory('D','D',size(fst),'gradient')
      deallocate(fst)
      call memory('D','D',size(bux1),'gradient')
      deallocate(bux1)
      call memory('D','D',size(bux2),'gradient')
      deallocate(bux2)
#ifdef MPI
      call memory('D','D',size(cl),'gradient')
      deallocate(cl)
      call memory('D','D',size(gl),'gradient')
      deallocate(gl)
      call memory('D','D',size(ftl2),'gradient')
      deallocate(ftl2)
      call memory('D','I',size(listftl),'gradient')
      deallocate(listftl)
#endif

      call timer('gradient',2)
      return
      end

