      subroutine denmat(c,eta,h,s,enum,
     .                 nbasis,nbands,ncmax,nctmax,
     .                 nfmax,nftmax,nhmax,
     .                 numc,listc,numct,listct,cttoc,
     .                 numf,listf,numft,listft,fttof,
     .                 numh,listh,
     .                 dm,edm,nbasisloc)
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
C real*8 h(nhmax,nbasisloc)   : Hamiltonian matrix (sparse)
C real*8 s(nhmax,nbasisloc)   : Overlap matrix (sparse)
C real*8 enum                 : Total number of electrons
C integer nbasis              : Number of atomic orbitals
C integer nbands              : Number of Localized Wave Functions
C integer ncmax               : First dimension of listc and C, and maximum
C                               number of nonzero elements of each row of C
C integer nctmax              : Max num of <>0 elements of each col of C
C integer nfmax               : Max num of <>0 elements of each row of 
C                               F = Ct x H
C integer nftmax              : Max num of <>0 elements of each col of F
C integer nhmax               : First dimension of listh and H, and maximum
C                               number of nonzero elements of each row of H
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
C integer numft(nbasis)        : Control vector of F transpose matrix
C                               (number of <>0  elements of each col of F)
C integer listft(nfmax,nbasisloc) : Control vector of F transpose matrix
C                               (list of <>0  elements of each col of F)
C integer fttof(nfmax,nbasisloc)  : Map from Ft to F indexing
C integer numh(nbasis)        : Control vector of H matrix
C                               (number of nonzero elements of each row of H)
C integer listh(nhmax,nbasis) : Control vector of H matrix
C                               (list of nonzero elements of each row of H)
C ************************* OUTPUT **********************************
C real*8 dm(nhmax,nbasisloc)  : Density Matrix
C real*8 edm(nhmax,nbasisloc) : Energy density matrix
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
     .  nbasis,nbands,ncmax,nctmax,nfmax,nhmax,nbasisloc

      integer
     .  cttoc(nctmax,*),listc(ncmax,nbasis),
     .  listct(nctmax,*),listf(nfmax,*),
     .  listh(nhmax,nbasis),
     .  numc(nbasis),numct(*),numf(*),
     .  numh(nbasis)
     
      integer
     .  nftmax,fttof(nftmax,nbasisloc),listft(nftmax,nbasisloc),
     .  numft(nbasis)

      double precision
     .  c(ncmax,nbasisloc),dm(nhmax,nbasisloc),edm(nhmax,nbasisloc),
     .  enum,eta,
     .  h(nhmax,nbasisloc),s(nhmax,nbasisloc)
     
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
     .  i, in, im, j, jm, jn, m, mi, mn, n, ni, nm,
     .  Node, Nodes, ii, nl, il, nn,
     .  nbandsmin, nbandsmax, nbandspernode, iloc, nbl,
     .  nbandsloc
#ifdef MPI
      integer
     .  MPIerror, nremainder, nloc, nbmin, nbmax, ml
      double precision
     .  rtmp
      double precision, dimension(:,:), allocatable, save :: 
     .  cl, hl, sl, chccl, csccl
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

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
      call GetNodeOrbs(nbasis,0,Nodes,nbl)
      allocate(cl(ncmax,nbl))
      call memory('A','D',ncmax*nbl,'denmat')
      allocate(hl(nhmax,nbl))
      call memory('A','D',nhmax*nbl,'denmat')
      allocate(sl(nhmax,nbl))
      call memory('A','D',nhmax*nbl,'denmat')

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

C Initialize temporary arrays ..............................................
      do nn = 1,nbandspernode+1
        do m = 1,nbasis
          cHrow(m,nn) = 0.d0
          cSrow(m,nn) = 0.d0
        enddo
      enddo
      do nn = 1,nbandspernode+1
        do i = 1,nbands
          cscrow(i,nn) = 0.d0
          chcrow(i,nn) = 0.d0
        enddo
      enddo

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

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do i = 1,nbasisloc
            do ii = 1,ncmax
              cl(ii,i) = c(ii,i)
            enddo
            do ii = 1,nhmax
              hl(ii,i) = h(ii,i)
              sl(ii,i) = s(ii,i)
            enddo
          enddo
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(hl(1,1),nhmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)
        call MPI_Bcast(sl(1,1),nhmax*nbl,DAT_double,nn-1,
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
#else
            im = cttoc(mi,i)
            cim = c(im,m)
#endif
            do nm = 1,numh(m)
              n = listh(nm,m)
#ifdef MPI
              Smn = sl(nm,ml)
              Hmn = hl(nm,ml)
#else
              Smn = S(nm,m)
              Hmn = H(nm,m)
#endif
              cSrow(n,i) = cSrow(n,i) + cim * Smn
              cHrow(n,i) = cHrow(n,i) + cim * Hmn
            enddo
#ifdef MPI
          endif
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
          do i = 1,nbasisloc
            do ii = 1,ncmax
              cl(ii,i) = c(ii,i)
            enddo
          enddo
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
          cscrow(i,iloc) = cscrow(i,iloc) + 2.d0
          chcrow(i,iloc) = chcrow(i,iloc) + 2.d0 * eta
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
          do i = 1,nbasisloc
            do ii = 1,ncmax
              cl(ii,i) = c(ii,i)
            enddo
          enddo
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
      allocate(chccl(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')
      allocate(csccl(nbasis,nbandspernode+1))
      call memory('A','D',nbasis*(nbandspernode+1),'denmat')
#endif

C Allocate column vectors and initialise
      allocate(chccCol(nbands,nbl))
      call memory('A','D',nbands*nbl,'denmat')
      allocate(csccCol(nbands,nbl))
      call memory('A','D',nbands*nbl,'denmat')

      do nl = 1,nbl
        do i = 1,nbands
          csccCol(i,nl) = 0.d0
          chccCol(i,nl) = 0.d0
        enddo
      enddo

C Loop over Nodes
      do nn = 1,Nodes

#ifdef MPI
C Broadcast chcc/cscc
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nloc = (nbmax-nbmin+1)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do i = 1,nloc
            do ii = 1,nbasis
              chccl(ii,i) = chcc(ii,i)
              csccl(ii,i) = cscc(ii,i)
            enddo
          enddo
        endif

C Broadcast arrays
        call MPI_Bcast(chccl(1,1),nbasis*nloc,DAT_double,
     .        nn-1,MPI_Comm_World,MPIerror)
        call MPI_Bcast(csccl(1,1),nbasis*nloc,DAT_double,
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
            if (i.ge.nbmin.and.i.le.nbmax) then
              iloc = i - nbmin + 1
              csccCol(i,n) = csccl(ni,iloc)
              chccCol(i,n) = chccl(ni,iloc)
            endif
#else
            csccCol(i,n) = cscc(ni,i)
            chccCol(i,n) = chcc(ni,i)
#endif
          enddo
        enddo

C End loop over Nodes
      enddo

#ifdef MPI
C Deallocate chccl and csccl as these are no longer needed
      call memory('D','D',size(chccl),'denmat')
      deallocate(chccl)
      call memory('D','D',size(csccl),'denmat')
      deallocate(csccl)
#endif

C Zero dm and edm
      do n = 1,nbasisloc
        nl = iptrltog(n)
        do mn = 1,numh(nl)
          dm(mn,n) = 0.0
          edm(mn,n) = 0.0
        enddo
      enddo

C Loop over Nodes
      do nn = 1,Nodes

#ifdef MPI
C Get local size on broadcast node
        call GetNodeOrbs(nbasis,nn-1,Nodes,nbl)

C Place arrays in local storage on broadcast node
        if (Node.eq.nn-1) then
          do i = 1,nbasisloc
            do ii = 1,ncmax
              cl(ii,i) = c(ii,i)
            enddo
          enddo
        endif

C Broadcast arrays
        call MPI_Bcast(cl(1,1),ncmax*nbl,DAT_double,nn-1,
     .        MPI_Comm_World,MPIerror)

C Build pointer to local orbitals for current working Node
        do m = 1,nbasis
          call GlobalToLocalOrb(m,nn-1,Nodes,ml)
          iptrgtol(m) = ml
        enddo
#endif

        do n = 1,nbasisloc
          nl = iptrltog(n)

C Find column n of c*cscc and c*chcc
C Use that H is symmetric to determine required elements
          do mn = 1,numh(nl)
            m = listh(mn,nl)
#ifdef MPI
            ml = iptrgtol(m)
            if (ml.gt.0) then
C Find element (m,n) of c*cscc and c*chcc
              ccsccmn = 0.d0
              cchccmn = 0.d0
              do im = 1,numc(m)
                i = listc(im,m)
                ccsccmn = ccsccmn + cl(im,ml) * csccCol(i,n)
                cchccmn = cchccmn + cl(im,ml) * chccCol(i,n)
              enddo
C Use fact that dm and edm are symmetric
              dm(mn,n)  = dm(mn,n) + 2.d0 * ccsccmn
              edm(mn,n) = edm(mn,n) + 2.d0 * cchccmn
            endif
#else
C Find element (m,n) of c*cscc and c*chcc
            ccsccmn = 0.d0
            cchccmn = 0.d0
            do im = 1,numc(m)
              i = listc(im,m)
              ccsccmn = ccsccmn + c(im,m) * csccCol(i,n)
              cchccmn = cchccmn + c(im,m) * chccCol(i,n)
            enddo
C Use fact that dm and edm are symmetric
            dm(mn,n)  = dm(mn,n) + 2.d0 * ccsccmn
            edm(mn,n) = edm(mn,n) + 2.d0 * cchccmn
#endif
          enddo
        enddo

C End loop over Nodes
      enddo

C ........................

C Normalize DM to exact charge .........................
C Calculate total output charge ...
      qout = 0.0d0
      do i = 1,nbasisloc
        ii = iptrltog(i)
        do in = 1,numh(ii)
          qout = qout + dm(in,i) * s(in,i)
        enddo
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
        do i = 1,nbasisloc
          ii = iptrltog(i)
          do in = 1,numh(ii)
            dm(in,i) = dm(in,i) * fact
            edm(in,i) = edm(in,i) * fact
          enddo
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
#endif

C Stop time counter and return ..................
      call timer('denmat',2)
      end

