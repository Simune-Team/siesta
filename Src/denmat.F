      subroutine denmat(c,eta,h,s,enum,nbasis,nbands,ncmax,nctmax,
     .                  nfmax,nftmax,maxnh,numc,listc,numct,listct,
     .                  cttoc,numf,listf,numft,listft,fttof,numh,
     .                  listhptr,listh,dm,edm,nbasisloc)
C *******************************************************************
C Subroutine to compute the Density and Energy Density matrices
C for the Order-N functional of Kim et al. (PRB 52, 1640 (95))
C (generalization of that proposed by Mauri et al, and Ordejon et al)
C
C Density Matrix:
C  D_mu,nu = 2 * C_i,mu * ( 2 * delta_i,j - S_i,j) * C_j,nu
C
C Energy Density Matrix:
C  E_mu,nu = 2 * C_i,mu * ( H_i,j + 2 * eta * (delta_i,j - S_i,j) ) * C_j,nu
C
C (The factor 2 is for spin)
C
C (See Ordejon et al, PRB 51, 1456 (95))
C
C The DM is normalized to the exact number of electrons!!!
C
C Written by P.Ordejon, Noviembre'96
C Modified by J.M.Soler, May'97
C Parallelisation introduced by J.D. Gale, April'99
C ************************** INPUT **********************************
C real*8 c(ncmax,nbasisloc)   : Localized Wave Functions (sparse)
C real*8 eta                  : Fermi level parameter of Kim et al.
C real*8 h(maxnh)             : Hamiltonian matrix (sparse)
C real*8 s(maxnh)             : Overlap matrix (sparse)
C real*8 enum                 : Total number of electrons
C integer nbasis              : Number of atomic orbitals
C integer nbands              : Number of Localized Wave Functions
C integer ncmax               : First dimension of listc and C, and maximum
C                               number of nonzero elements of each row of C
C integer nctmax              : Max num of <>0 elements of each col of C
C integer nfmax               : Max num of <>0 elements of each row of 
C                               F = Ct x H
C integer nftmax              : Max num of <>0 elements of each col of F
C integer maxnh               : First dimension of listh and H, and maximum
C                               number of nonzero elements of H
C integer numc(nbasis)        : Control vector of C matrix
C                               (number of nonzero elements of each row of C)
C integer listc(ncmax,nbasis) : Control vector of C matrix
C                              (list of nonzero elements of each row of C)
C integer numct(nbands)       : Control vector of C transpose matrix
C                              (number of <>0  elements of each col of C)
C integer listct(ncmax,nbands): Control vector of C transpose matrix
C                              (list of <>0  elements of each col of C)
C integer cttoc(ncmax,nbands) : Map from Ct to C indexing
C integer numf(nbands)        : Control vector of F matrix
C                               (number of <>0  elements of each row of F)
C integer listf(nfmax,nbands) : Control vector of F matrix
C                               (list of <>0  elements of each row of F)
C integer numft(nbasis)       : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,nbasisloc) : Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer fttof(nfmax,nbasisloc)  : Map from Ft to F indexing
C integer numh(nbasisloc)     : Control vector of H matrix
C                               (number of nonzero elements of each row of H)
C integer listhptr(nbasisloc) : Control vector of H matrix
C                               (pointer to the start of each row of H)
C integer listh(maxnh)        : Control vector of H matrix
C                               (list of nonzero elements of each row of H)
C ************************* OUTPUT **********************************
C real*8 dm(maxnh)            : Density Matrix
C real*8 edm(maxnh)           : Energy density matrix
C *******************************************************************
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
     .  nbasis,nbands,ncmax,nctmax,nfmax,maxnh,nbasisloc

      integer
     .  cttoc(nctmax,*),listc(ncmax,nbasis),listct(nctmax,*),
     .  listf(nfmax,*),listh(maxnh),numc(nbasis),numct(*),
     .  numf(*),numh(nbasisloc),listhptr(nbasisloc)
     
      integer
     .  nftmax,fttof(nftmax,nbasisloc),listft(nftmax,nbasisloc),
     .  numft(nbasis)

      double precision
     .  c(ncmax,nbasisloc),dm(maxnh),edm(maxnh),enum,eta,
     .  h(maxnh),s(maxnh)
     
      external
     .  timer, memory

C Internal variables .................................................
C   Notation hints:
C     m,n : basis orbital inexes (mu,nu)
C     i,j : band (and LWF) indexes
C     im  : index for LWF's of basis orbital m
C     mi  : index for basis orbitals of LWF i
C     nm  : index for basis orbitals connected to basis orbital m

      integer 
     .  i, ii, in, ind, im, il, iloc, j, jm, jn, m, mi, mn, n, nh, ni,
     .  Node, Nodes, nl, nm, nn, nbandsmin, nbandsmax, nbandspernode,
     .  nbl, nbandsloc
#ifdef MPI
      integer
     .  MPIerror, ml, nremainder, nloc, nbmin, nbmax, maxnhg
      integer, dimension(:), allocatable ::
     .  listhlptr, listhl, numhl
      double precision
     .  rtmp
      double precision, dimension(:), allocatable, save :: 
     .  hl, sl, chccl
      double precision, dimension(:,:), allocatable, save :: 
     .  cl
#endif
      integer, dimension(:), allocatable, save :: 
     .  iptrgtol, iptrltog

      double precision, dimension(:,:), allocatable, save :: 
     .  cHrow, cSrow, chcrow, cscrow, chccCol, csccCol,
     .  chcc, cscc

      double precision
     .  cim, cnj, chin, csin, chccim, csccim, cchccmn, ccsccmn,
     .  Hmn, Smn, qout, fact, cSn, cHn

C ........................

C Start time counter .....................................................
      call timer('denmat',1)

#ifdef MPI
C Get Node number
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)

C Get global sizes for first dimension of H/S
      call MPI_AllReduce(maxnh,maxnhg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Allocate local memory for c/h/s
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(cl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'denmat')
      allocate(hl(maxnhg))
      call memory('A','D',maxnhg,'denmat')
      allocate(sl(maxnhg))
      call memory('A','D',maxnhg,'denmat')

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

C .......................

C Allocate workspace arrays
      allocate(cHrow(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')
      allocate(cSrow(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')
      allocate(chcrow(nbands,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')
      allocate(cscrow(nbands,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')

#ifdef MPI
C Allocate local storage for H/S control vectors
      call GetNodeOrbs(nbasis,Node,Nodes,nbl)
      allocate(numhl(nbasis))
      call memory('A','I',nbasis,'denmat')
      allocate(listhlptr(nbasis))
      call memory('A','I',nbasis,'denmat')
      allocate(listhl(maxnhg))
      call memory('A','I',maxnhg,'denmat')
#endif

C Initialize temporary arrays ..............................................
      cHrow(1:nbasis,1:nbandspernode+1) = 0.d0
      cSrow(1:nbasis,1:nbandspernode+1) = 0.d0
      cscrow(1:nbands,1:nbandspernode+1) = 0.d0
      chcrow(1:nbands,1:nbandspernode+1) = 0.d0

C Build pointer for local to global
      allocate(iptrgtol(nbasis))
      call memory('A','I',nbasis,'denmat')
      allocate(iptrltog(nbasis))
      call memory('A','I',nbasis,'denmat')
      do i = 1,nbasisloc
        call LocalToGlobalOrb(i,Node,Nodes,il)
        iptrltog(i) = il
      enddo

C ........................

C Loop over Nodes for broadcasting of c,h and s
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
#endif

C Find cscc=(2-ct*S*c)*ct and chcc=(ct*H*c+2eta(1-ct*S*c))*ct.
        do i = 1,nbandsloc
      
C       Find row i of cS=ct*S and cH=ct*H
          do mi = 1,numct(i)
            m = listct(mi,i)
#ifdef MPI
            call GlobalToLocalOrb(m,nn-1,Nodes,ml)
            if (ml.gt.0) then
              im = cttoc(mi,i)
              cim = cl(im,ml)
              do nm = 1,numhl(ml)
                ind = listhlptr(ml)+nm
                n = listhl(ind)
                Smn = sl(ind)
                Hmn = hl(ind)
                cSrow(n,i) = cSrow(n,i) + cim * Smn
                cHrow(n,i) = cHrow(n,i) + cim * Hmn
              enddo
            endif
#else
            im = cttoc(mi,i)
            cim = c(im,m)
            do nm = 1,numh(m)
              ind = listhptr(m)+nm
              n = listh(ind)
              Smn = S(ind)
              Hmn = H(ind)
              cSrow(n,i) = cSrow(n,i) + cim * Smn
              cHrow(n,i) = cHrow(n,i) + cim * Hmn
            enddo
#endif
          enddo

        enddo

C End loop over Nodes
      enddo

#ifdef MPI
C Deallocate hl and sl as these arrays are no longer needed
      call memory('D','D',size(hl),'denmat')
      deallocate(hl)
      call memory('D','D',size(sl),'denmat')
      deallocate(sl)
#endif

C Loop over Nodes for broadcasting of c
      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          cl(1:ncmax,1:nbasisloc) = c(1:ncmax,1:nbasisloc)
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

      do i = nbandsmin,nbandsmax
        iloc = i - nbandsmin + 1

C Find row i of csc=2-ct*S*c and chc=ct*H*c+2eta(1-ct*S*c)
C Now use the list of nonzero elements of f=ct*H
        do ni = 1,numf(iloc)
          n = listf(ni,iloc)
#ifdef MPI
          call GlobalToLocalOrb(n,nn-1,Nodes,nl)
          if (nl.gt.0) then
#endif
            cSn = cSrow(n,iloc)
            cHn = cHrow(n,iloc)
            csin = - cSn
            chin = cHn - 2.d0*eta*cSn
            do jn = 1,numc(n)
              j = listc(jn,n)
#ifdef MPI
              cnj = cl(jn,nl)
#else
              cnj = c(jn,n)
#endif
              cscrow(j,iloc) = cscrow(j,iloc) + csin * cnj
              chcrow(j,iloc) = chcrow(j,iloc) + chin * cnj
            enddo
#ifdef MPI
          endif
#endif

        enddo

C       Add on diagonal terms 2 and 2eta
        if (Node.eq.nn-1) then
          cscrow(i,iloc) = cscrow(i,iloc) + 2.0d0
          chcrow(i,iloc) = chcrow(i,iloc) + 2.0d0 * eta
        endif

      enddo

C End loop over Nodes
      enddo

C Deallocate cHrow/cSrow as they are no longer needed
      call memory('D','D',size(cHrow),'denmat')
      deallocate(cHrow)
      call memory('D','D',size(cSrow),'denmat')
      deallocate(cSrow)

C Allocate chcc/cscc
      allocate(chcc(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')
      allocate(cscc(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')
        
C Loop over Nodes for broadcasting of c
      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          cl(1:ncmax,1:nbasisloc) = c(1:ncmax,1:nbasisloc)
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
#endif

      do i = nbandsmin,nbandsmax
        iloc = i - nbandsmin + 1

C       Find row i of cscc=csc*ct and chcc=chc*ct. 
C       Only the nonzero elements of f=cH will be required.
        do mi = 1,numf(iloc)
          m = listf(mi,iloc)
#ifdef MPI
          call GlobalToLocalOrb(m,nn-1,Nodes,ml)
          if (ml.gt.0) then
#endif
            csccim = 0.d0
            chccim = 0.d0
            do jm = 1,numc(m)
              j = listc(jm,m)
#ifdef MPI
              csccim = csccim + cscrow(j,iloc) * cl(jm,ml)
              chccim = chccim + chcrow(j,iloc) * cl(jm,ml)
#else
              csccim = csccim + cscrow(j,iloc) * c(jm,m)
              chccim = chccim + chcrow(j,iloc) * c(jm,m)
#endif
            enddo
#ifdef MPI
            cscc(mi,iloc) = csccim
            chcc(mi,iloc) = chccim
          endif
#else
          cscc(mi,i) = csccim
          chcc(mi,i) = chccim
#endif
        enddo

      enddo

C End loop over Nodes
      enddo

C Deallocate chcrow/cscrow as they are no longer needed
      call memory('D','D',size(chcrow),'denmat')
      deallocate(chcrow)
      call memory('D','D',size(cscrow),'denmat')
      deallocate(cscrow)

C ........................

      call GetNodeOrbs(nbasis,0,Nodes,nbl)
#ifdef MPI
C Allocate workspace for dm and edm on local node
      allocate(chccl(nbasis))
      call memory('A','D',nbasis,'denmat')
#endif

C Allocate column vectors and initialise
      allocate(chccCol(nbands,nbl))
      call memory('A','D',nbands*nbl,'denmat')
      allocate(csccCol(nbands,nbl))
      call memory('A','D',nbands*nbl,'denmat')

      csccCol(1:nbands,1:nbl) = 0.0d0
      chccCol(1:nbands,1:nbl) = 0.0d0

C Loop over Nodes
      do nn = 1,Nodes

#ifdef MPI
C Broadcast chcc
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nloc = (nbmax-nbmin+1)

        do iloc = 1,nloc
C Place arrays in local storage on broadcast node
          if (Node.eq.nn-1) then
            chccl(1:nbasis) = chcc(1:nbasis,iloc)
          endif

C Broadcast arrays
          call MPI_Bcast(chccl,nbasis,DAT_double,
     .        nn-1,MPI_Comm_World,MPIerror)
#endif

C Find dm=c*cscc and edm=c*chcc. Only the nonzero elements of H.
          do n = 1,nbasisloc
            nl = iptrltog(n)

C Use listft to expand a column of cscc
            do in = 1,numft(nl)
              i = listft(in,n) 
              ni = fttof(in,n)
#ifdef MPI
              if (i.eq.(iloc+nbmin-1)) then
                chccCol(i,n) = chccl(ni)
              endif
#else
              csccCol(i,n) = cscc(ni,i)
              chccCol(i,n) = chcc(ni,i)
#endif
            enddo
          enddo
#ifdef MPI
        enddo
#endif

C End loop over Nodes
      enddo

#ifdef MPI
C Loop over Nodes
      do nn = 1,Nodes

C Broadcast cscc
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nloc = (nbmax-nbmin+1)

        do iloc = 1,nloc
C Place arrays in local storage on broadcast node
          if (Node.eq.nn-1) then
            chccl(1:nbasis) = cscc(1:nbasis,iloc)
          endif

C Broadcast arrays
          call MPI_Bcast(chccl,nbasis,DAT_double,
     .        nn-1,MPI_Comm_World,MPIerror)

C Find dm=c*cscc and edm=c*chcc. Only the nonzero elements of H.
          do n = 1,nbasisloc
            nl = iptrltog(n)

C Use listft to expand a column of cscc
            do in = 1,numft(nl)
              i = listft(in,n)
              ni = fttof(in,n)
              if (i.eq.(iloc+nbmin-1)) then
                csccCol(i,n) = chccl(ni)
              endif
            enddo
          enddo
        enddo

C End loop over Nodes
      enddo

C Deallocate chccl as this is no longer needed
      call memory('D','D',size(chccl),'denmat')
      deallocate(chccl)
#endif

C Zero dm and edm
      nh = listhptr(nbasisloc) + numh(nbasisloc)
      dm(1:nh) = 0.0d0
      edm(1:nh) = 0.0d0

C Loop over Nodes
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
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .    MPI_Comm_World,MPIerror)

C Build pointer to local orbitals for current working Node
        do m = 1,nbasis
          call GlobalToLocalOrb(m,nn-1,Nodes,ml)
          iptrgtol(m) = ml
        enddo
#endif

        do n = 1,nbasisloc

C Find column n of c*cscc and c*chcc
C Use that H is symmetric to determine required elements
          do mn = 1,numh(n)
            ind = listhptr(n)+mn
            m = listh(ind)
#ifdef MPI
            ml = iptrgtol(m)
            if (ml.gt.0) then
C Find element (m,n) of c*cscc and c*chcc
              ccsccmn = 0.0d0
              cchccmn = 0.0d0
              do im = 1,numc(m)
                i = listc(im,m)
                ccsccmn = ccsccmn + cl(im,ml) * csccCol(i,n)
                cchccmn = cchccmn + cl(im,ml) * chccCol(i,n)
              enddo
C Use fact that dm and edm are symmetric
              dm(ind)  = dm(ind) + 2.0d0 * ccsccmn
              edm(ind) = edm(ind) + 2.0d0 * cchccmn
            endif
#else
C Find element (m,n) of c*cscc and c*chcc
            ccsccmn = 0.0d0
            cchccmn = 0.0d0
            do im = 1,numc(m)
              i = listc(im,m)
              ccsccmn = ccsccmn + c(im,m) * csccCol(i,n)
              cchccmn = cchccmn + c(im,m) * chccCol(i,n)
            enddo
C Use fact that dm and edm are symmetric
            dm(ind)  = dm(ind) + 2.0d0 * ccsccmn
            edm(ind) = edm(ind) + 2.0d0 * cchccmn
#endif
          enddo
        enddo

C End loop over Nodes
      enddo

C ........................

C Normalize DM to exact charge .........................
C Calculate total output charge ...
      qout = 0.0d0
      nh = listhptr(nbasisloc) + numh(nbasisloc)
      do i = 1,nh
        qout = qout + dm(i) * s(i)
      enddo
#ifdef MPI
      call MPI_AllReduce(qout,rtmp,1,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      qout = rtmp
#endif
      if (Node.eq.0) then
      write(6,"(/a,f12.4)") 'denmat: qtot (before DM normalization) = ',
     .              qout
      endif
C ...

      if (dabs(enum-qout) .gt. 0.05d0) then
        fact = enum / qout
      
C Normalize ...
        do i = 1,nh
          dm(i) = dm(i) * fact
          edm(i) = edm(i) * fact
        enddo
C ...
      endif
C ........................

C Free local memory
      call memory('D','I',size(iptrgtol),'denmat')
      deallocate(iptrgtol)
      call memory('D','I',size(iptrltog),'denmat')
      deallocate(iptrltog)
      call memory('D','D',size(chcc),'denmat')
      deallocate(chcc)
      call memory('D','D',size(cscc),'denmat')
      deallocate(cscc)
      call memory('D','D',size(chccCol),'denmat')
      deallocate(chccCol)
      call memory('D','D',size(csccCol),'denmat')
      deallocate(csccCol)
#ifdef MPI
      call memory('D','D',size(cl),'denmat')
      deallocate(cl)
      call memory('D','I',size(numhl),'denmat')
      deallocate(numhl)
      call memory('D','I',size(listhlptr),'denmat')
      deallocate(listhlptr)
      call memory('D','I',size(listhl),'denmat')
      deallocate(listhl)
#endif

C Stop time counter and return ..................
      call timer('denmat',2)
      end

