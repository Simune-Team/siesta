      subroutine gradient(c,eta,enum,h,s,nbasis,nbands,ncmax,nctmax,
     .                    nfmax,nftmax,maxnh,nhijmax,numc,listc,
     .                    numct,listct,cttoc,numf,listf,numft,listft,
     .                    numh,listhptr,listh,numhij,listhij,grad,
     .                    ener,nbasisloc)

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
C integer listft(nfmax,nbasisloc) : Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer numh(nbasisloc)      : Control vector of H matrix
C                                (number of <>0  elements of each row of H)
C integer listhptr(nbasisloc)  : Control vector of H matrix
C                                (pointer to the start of each row of H)
C integer listh(maxnh)         : Control vector of H matrix 
C                               (list of <>0  elements of each row of H)
C integer numhij(nbasis)       : Control vector of Hij matrix
C                                (number of <>0  elements of each row of Hij)
C integer listhij(nhijmax,nbands): Control vector of Hij matrix 
C                                (list of <>0  elements of each row of Hij)
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
      use mpi
#endif

      implicit none

      integer
     .  nbasis,nbands,ncmax,nctmax,nfmax,nftmax,maxnh,nhijmax,nbasisloc

      integer
     .  cttoc(nctmax,*),listc(ncmax,nbasis),listct(nctmax,*),
     .  listf(nfmax,*),listft(nftmax,nbasisloc),listh(maxnh),
     .  listhij(nhijmax,*),numc(nbasis),numct(*),numf(*),
     .  numft(nbasis),numh(nbasisloc),listhptr(nbasisloc),
     .  numhij(*)

      double precision
     .  c(ncmax,nbasisloc),ener,eta,enum,grad(ncmax,nbasisloc),
     .  h(maxnh),s(maxnh)

C Internal variables ......................................................

      integer
     .  i,ik,in,ind,j,jn,k,kl,kn,knlh,mu,muk,Node,Nodes,
     .  nn, nbandsmin,nbandsmax,nbandspernode, iloc, 
     .  nbandsloc
      integer, dimension(:), allocatable, save :: iptr
#ifdef MPI
      integer
     .  MPIerror, ii, il, jloc, mul, nh, nremainder, 
     .  nbmin, nbmax, maxnhg, nloc, nbl
      integer, dimension(:), allocatable ::
     .  listhlptr, listhl, numhl
      integer, dimension(:,:), allocatable, save :: 
     .  listftl
      double precision
     .  rtmp1(2),rtmp2(2)
      double precision, dimension(:), allocatable, save :: 
     .  hl, sl, ftl
      double precision, dimension(:,:), allocatable, save :: 
     .  cl, gl, ftl2
#endif
      double precision, dimension(:,:), allocatable, save :: 
     .  f, ft, fs, fst, bux1, bux2

      double precision
     .  a0,b0,p0,func1,func2
     
C..................

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
#else
      Node = 0
      Nodes = 1
      nbandsmin = 1
      nbandsmax = nbands
      nbandspernode = nbands
      nbandsloc = nbands
#endif

      call timer('gradient',1)
      
C Allocate dynamic memory
      allocate(iptr(nbasis))
      call memory('A','I',nbasis,'gradient')

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

      allocate(f(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'gradient')
      allocate(fs(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'gradient')
      allocate(bux1(nbands,nbandspernode+1))
      call memory('A','D',nbands*(nbandspernode+1),'gradient')
      allocate(bux2(nbands,nbandspernode+1))
      call memory('A','D',nbands*(nbandspernode+1),'gradient')

      f(1:nbasis,1:nbandspernode+1) = 0.0d0
      fs(1:nbasis,1:nbandspernode+1) = 0.0d0
      bux1(1:nbands,1:nbandspernode+1) = 0.0d0
      bux2(1:nbands,1:nbandspernode+1) = 0.0d0

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
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(hl,maxnhg,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(sl,maxnhg,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)

C Build pointer to local orbitals for current working Node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
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
              do kn = 1,numhl(kl)
                ind = listhlptr(kl)+kn
                knlh = listhl(ind)
                f(knlh,i) = f(knlh,i) + p0 * (hl(ind)-eta*sl(ind))
                fs(knlh,i) = fs(knlh,i) + p0 * sl(ind)
              enddo
            endif
#else
            ik = cttoc(in,i)
            p0 = c(ik,k)
            do kn = 1,numh(k)
              ind = listhptr(k)+kn
              knlh = listh(ind)
              f(knlh,i) = f(knlh,i) + p0 * ( h(ind) - eta*s(ind) )
              fs(knlh,i) = fs(knlh,i) + p0 * s(ind)
            enddo
#endif
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
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)

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

          do in = 1,numf(i)
            k = listf(in,i)
#ifdef MPI
            kl = iptr(k)
            if (kl.gt.0) then
#endif
              a0 = f(k,i)
              b0 = fs(k,i)
              do kn = 1,numc(k)
#ifdef MPI
                bux1(listc(kn,k),i) = bux1(listc(kn,k),i) + 
     .            a0 * cl(kn,kl)
                bux2(listc(kn,k),i) = bux2(listc(kn,k),i) + 
     .            b0 * cl(kn,kl)
#else
                bux1(listc(kn,k),i) = bux1(listc(kn,k),i) + 
     .            a0 * c(kn,k)
                bux2(listc(kn,k),i) = bux2(listc(kn,k),i) + 
     .            b0 * c(kn,k)
#endif
              enddo
#ifdef MPI
            endif
#endif
          enddo

        enddo

C End loop over Nodes
      enddo

      do i = nbandsmin,nbandsmax
        iloc = i - nbandsmin + 1
        func1 = func1 + bux1(i,iloc)
        do jn = 1,numhij(iloc)
          j = listhij(jn,iloc)
          func2 = func2 + bux1(j,iloc) * bux2(j,iloc)
        enddo
      enddo

C Allocate memory for transpose array fst
      allocate(fst(nbands,nbasisloc))
      call memory('A','D',nbands*nbasisloc,'gradient')
#ifdef MPI
      allocate(ftl(nbasis))
      call memory('A','D',nbasis,'gradient')
#endif

C Build pointer to local orbitals for current working Node
      do k = 1,nbasis
        call GlobalToLocalOrb(k,Node,Nodes,kl)
        iptr(k) = kl
      enddo

C Calculate transpose of fs
      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nloc = (nbmax-nbmin+1)

        do i = 1,nloc
C Place arrays in local storage on broadcast node
          if (Node.eq.nn-1) then
            do ii = 1,nbasis
              ftl(ii) = fs(ii,i)
            enddo
          endif

C Broadcast arrays
          call MPI_Bcast(ftl,nbasis,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)

          iloc = i + nbmin - 1
          do k = 1,nbasis
            kl = iptr(k)
            if (kl.gt.0) then
              fst(iloc,kl) = ftl(k)
            endif
          enddo
        enddo
#else
        do i = nbandsmin,nbandsmax
          iloc = i - nbandsmin + 1
          do k = 1,nbasis
            kl = iptr(k)
            if (kl.gt.0) then
              fst(i,k) = fs(k,iloc)
            endif
          enddo
        enddo
#endif
      enddo

C Deallocate fs as it is not used again
      call memory('D','D',size(fs),'gradient')
      deallocate(fs)

C Allocate memory for transpose array ft
      allocate(ft(nbands,nbasisloc))
      call memory('A','D',nbands*nbasisloc,'gradient')

C Calculate transpose of f
      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nloc = (nbmax-nbmin+1)

        do i = 1,nloc
C Place arrays in local storage on broadcast node
          if (Node.eq.nn-1) then
            do ii = 1,nbasis
              ftl(ii) = f(ii,i)
            enddo
          endif

C Broadcast arrays
          call MPI_Bcast(ftl,nbasis,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)

          iloc = i + nbmin - 1
          do k = 1,nbasis
            kl = iptr(k)
            if (kl.gt.0) then
              ft(iloc,kl) = ftl(k)
            endif
          enddo
        enddo
#else
        do i = nbandsmin,nbandsmax
          iloc = i - nbandsmin + 1
          do k = 1,nbasis
            kl = iptr(k)
            if (kl.gt.0) then
              ft(i,k) = f(k,iloc)
            endif
          enddo
        enddo
#endif
      enddo

#ifdef MPI
C Free ftl with old dimensions 
      call memory('D','D',size(ftl),'gradient')
      deallocate(ftl)

C Add contribution from f to grad
      do nn = 1,Nodes

C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)
        cl(1:ncmax,1:nbl) = 0.0

        do kl = 1,nbl
          call LocalToGlobalOrb(kl,nn-1,Nodes,k)
          do ik = 1,numc(k)
            j = listc(ik,k)
            if (j.ge.nbandsmin.and.j.le.nbandsmax) then
              jloc = j - nbandsmin + 1
              cl(ik,kl) = 4.d0 * f(k,jloc) + cl(ik,kl)
            endif
          enddo
        enddo

C Global reduction of grad
        call MPI_Reduce(cl(1,1),grad(1,1),ncmax*nbl,
     .    DAT_double,MPI_sum,nn-1,
     .    MPI_Comm_World,MPIerror)

      enddo
#else
      do k = 1,nbasis
        do ik = 1,numc(k)
          j = listc(ik,k)
          grad(ik,k) = 4.d0 * f(k,j) + grad(ik,k)
        enddo
      enddo
#endif

C Deallocate f
      call memory('D','D',size(f),'gradient')
      deallocate(f)

#ifdef MPI
C Create local storage memory for transpose arrays
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(ftl2(nbands,nbl))
      call memory('A','D',nbands*nbl,'gradient')
      allocate(listftl(nftmax,nbl))
      call memory('A','I',nftmax*nbl,'gradient')
      allocate(gl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'gradient')
#endif

      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)
        cl(1:ncmax,1:nbl) = 0.0

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          ftl2(1:nbands,1:nbasisloc) = ft(1:nbands,1:nbasisloc)
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
        call MPI_Bcast(ftl2(1,1),nbands*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(listftl(1,1),nftmax*nbl,MPI_integer,
     .        nn-1,MPI_Comm_World,MPIerror)

C Build pointer to local orbitals for current working Node
        do k = 1,nbasis
          call GlobalToLocalOrb(k,nn-1,Nodes,kl)
          iptr(k) = kl
        enddo
#endif

        do i = 1,nbandsloc

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
                a0 = a0 + bux2(j,i) * ftl2(j,mul)
              enddo
              cl(cttoc(ik,i),mul) = cl(cttoc(ik,i),mul) -2.0 * a0
            endif
          enddo
#else
          do ik = 1,numct(i)
            mu = listct(ik,i)
            a0 = 0.d0
            do muk = 1,numft(mu)
              j = listft(muk,mu)
              a0 = a0 + bux1(j,i) * fst(j,mu) + bux2(j,i) * 
     .          ft(j,mu)
            enddo
            grad(cttoc(ik,i),mu) = (-2.d0) * a0
          enddo
#endif
        enddo

#ifdef MPI
C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          ftl2(1:nbands,1:nbasisloc) = fst(1:nbands,1:nbasisloc)
        endif

C Broadcast arrays
        call MPI_Bcast(ftl2(1,1),nbands*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)

        do i = 1,nbandsloc

C multiply Hij x Fs and Sij x F row by row
C (only products of neccesary elements)
          do ik = 1,numct(i)
            mu = listct(ik,i)
            mul = iptr(mu)
            if (mul.gt.0) then
              a0 = 0.d0
              do muk = 1,numft(mu)
                j = listftl(muk,mul)
                a0 = a0 + bux1(j,i) * ftl2(j,mul)
              enddo
              cl(cttoc(ik,i),mul) = cl(cttoc(ik,i),mul) -2.0 * a0
            endif
          enddo
        enddo
#endif

#ifdef MPI
C Global reduction of grad
        call MPI_Reduce(cl(1,1),gl(1,1),ncmax*nbl,
     .    DAT_double,MPI_sum,nn-1,
     .    MPI_Comm_World,MPIerror)
        if (Node.eq.nn-1) then
          do k = 1,nbasisloc
            call LocalToGlobalOrb(kl,nn-1,Nodes,k)
            do ik = 1,numc(k)
              grad(ik,k) = gl(ik,k) + grad(ik,k)
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
      call MPI_AllReduce(rtmp1,rtmp2,2,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      func1 = rtmp2(1)
      func2 = rtmp2(2)
#endif

      ener = 2.d0 * func1 - func2 + eta * enum / 2.d0
C ...................

C Free dynamic memory
      call memory('D','I',size(iptr),'gradient')
      deallocate(iptr)
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

