      subroutine ctrans1(nr,ncmin,ncmax,nmax,ntmax,
     .              num,list,numt,listt,cttoc)
C ********************************************************************
C Finds the C transpose matrix control vectors numt and listt,
C and the index vector cttoc.
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nr             : Number of rows of C 
C integer ncmin          : First column of full C stored locally
C integer ncmax          : Last column of full C stored locally
C integer nc             : Right-hand dimension of output arrays
C integer nmax           : First dimension of list and C, and maximum
C                           number of nonzero elements of each row of C
C integer ntmax          : maximum number of nonzero elements of each
C                           column of C
C integer num(nr)        : Control vector of C matrix
C                           (number of nonzero elements of each row of C)
C integer list(nmax,nr)  : Control vector of C matrix
C                          (list of nonzero elements of each row of C)
C **************************** OUTPUT *********************************
C integer numt(ncloc)     : Control vector of C matrix
C                         (number of nonzero elements of each column of C)
C integer listt(ntmax,ncloc): Control vector of C transpose matrix
C                           (list of nonzero elements of each column of C)
C integer cttoc(ntmax,ncloc): Map from C transpose to C indexing
C *********************************************************************
      implicit none

      integer
     .  ncmin,ncmax,nmax,nr,ntmax,
     .  cttoc(ntmax,ncmax-ncmin+1),list(nmax,nr),
     .  listt(ntmax,ncmax-ncmin+1),
     .  num(nr),numt(ncmax-ncmin+1)

C Internal variables ..................................................
      integer
     .  i,imu,mu,n,il,ncloc
C ..........................
C  Initialize numt list ...............................................
      ncloc = ncmax - ncmin + 1
      do i = 1,ncloc
        numt(i) = 0
      enddo
C ..........................
C  Construct information for transpose of C ............................
      do mu = 1,nr
        do imu = 1,num(mu)
          i = list(imu,mu)
          if (i.ge.ncmin.and.i.le.ncmax) then
            il = i - ncmin + 1
            numt(il) = numt(il)+1
            n = numt(il)
            if (n .le. ntmax) then
              listt(n,il) = mu
              cttoc(n,il) = imu
            endif
          endif
        enddo
      enddo
C ..........................
C Check dimensions .....................................................
      do i = 1,ncloc
        call chkdim('ctrans','ntmax',ntmax,numt(i),1)
      enddo
C ..........................
      return
      end
      subroutine ctrans2(nr,nc,nmax,ntmax,
     .              num,list,numt,listt,cttoc,Node,Nodes)
C ********************************************************************
C Finds the C transpose matrix control vectors numt and listt,
C and the index vector cttoc.
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nr             : Number of rows of C 
C integer nc             : Number of columns of full C
C integer nmax           : First dimension of list and C, and maximum
C                           number of nonzero elements of each row of C
C integer ntmax          : maximum number of nonzero elements of each
C                           column of C
C integer num(nr)        : Control vector of C matrix
C                           (number of nonzero elements of each row of C)
C integer list(nmax,nr)  : Control vector of C matrix
C                          (list of nonzero elements of each row of C)
C integer Node           : Local node number
C integer Nodes          : Total number of global nodes
C **************************** OUTPUT *********************************
C integer numt(nc)       : Control vector of C matrix
C                         (number of nonzero elements of each column of C)
C integer listt(ntmax,nc) : Control vector of C transpose matrix
C                           (list of nonzero elements of each column of C)
C integer cttoc(ntmax,nc) : Map from C transpose to C indexing
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

      integer
     .  nc,nmax,nr,ntmax,Node,Nodes,
     .  cttoc(ntmax,*),list(nmax,nr),listt(ntmax,*),
     .  num(nr),numt(nc)

#ifdef MPI
      integer
     .  MPIerror, nbandspernode, nbandsmin, nbandsmax, nremainder,
     .  in,nrmax,nrmin,nrloc,mul
      integer, dimension(:,:), allocatable, save :: listlocal
      integer, dimension(:)  , allocatable, save :: numlocal
#endif

C Internal variables ..................................................
      integer
     .  i,imu,mu,n,il
C ..........................
C  Initialize numt list ...............................................
      do i = 1,nc
        numt(i) = 0
      enddo
C ..........................
C  Construct information for transpose of C ............................
C  Need to use information from other Nodes
#ifdef MPI
C
C  Work out local distribution of bands
C
      nbandspernode = ((nr-1)/Nodes)
      nremainder = nr - nbandspernode*Nodes
      nbandsmin = Node*nbandspernode + min(nremainder,Node) + 1
      nbandsmax = nbandsmin + nbandspernode - 1
      if (Node.lt.nremainder) nbandsmax = nbandsmax + 1
C
C  Allocate local memory
C
      allocate(numlocal(nbandspernode+1))
      call memory('A','I',nbandspernode+1,'ctrans2')
      allocate(listlocal(nmax,nbandspernode+1))
      call memory('A','I',nmax*(nbandspernode+1),'ctrans2')
C
C  Loop over Nodes broadcasting num/list
C
      do in = 0,Nodes-1
C
C  Work out which bands are on Node = in
C
        nrmin = in*nbandspernode + min(nremainder,in) + 1
        nrmax = nrmin + nbandspernode - 1
        if (in.lt.nremainder) nrmax = nrmax + 1
        nrloc = nrmax - nrmin + 1
C
C  Copy data on broadcast node to workspace arrays
C
        if (Node.eq.in) then
          do i = 1,nrloc
            numlocal(i) = num(i)
            do il = 1,nmax
              listlocal(il,i) = list(il,i)
            enddo
          enddo
        endif
C
C  Transfer data
C
        call MPI_Bcast(numlocal,nrloc,MPI_integer,in,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(listlocal(1,1),nmax*nrloc,MPI_integer,in,
     .    MPI_Comm_World,MPIerror)
C
C  Use information from Node = in
C
        do mu = nrmin,nrmax
          mul = mu - nrmin + 1
          do imu = 1,numlocal(mul)
            i = listlocal(imu,mul)
            numt(i) = numt(i)+1
            n = numt(i)
            if (n .le. ntmax) then
              call GlobalToLocalOrb(i,Node,Nodes,il)
              if (il.gt.0) then
                listt(n,il) = mu
                cttoc(n,il) = imu
              endif
            endif
          enddo
        enddo
      enddo
C
C  Free local memory
C
      call memory('D','I',size(listlocal),'ctrans2')
      deallocate(listlocal)
      call memory('D','I',size(numlocal),'ctrans2')
      deallocate(numlocal)
#else
      do mu = 1,nr
        do imu = 1,num(mu)
          i = list(imu,mu)
          numt(i) = numt(i)+1
          n = numt(i)
          if (n .le. ntmax) then
            call GlobalToLocalOrb(i,Node,Nodes,il)
            if (il.gt.0) then
              listt(n,il) = mu
              cttoc(n,il) = imu
            endif
          endif
        enddo
      enddo
#endif
C ..........................
C Check dimensions .....................................................
      do i = 1,nc
        call chkdim('ctrans','ntmax',ntmax,numt(i),1)
      enddo
C ..........................
      return
      end
      subroutine axb_build1(nramin,nramax,nca,namax,numa,lista,
     .               nrb,ncb,nbmax,numb,listbptr,listb,ncmax,
     .               numc,listc,node,Nodes)
C ********************************************************************
C Constructs control indexes of a C matrix in sparse form,
C C being the product of A and B (also in sparse form)
C
C              C = A x B
C
C In full form: A is rectangular, and has dimension:  nra x nca
C               B is rectangular, and has dimension:  nrb x ncb
C and, as a result:
C               C is rectangular, and has dimension:  nra x ncb
C (Of course, nca must be equal to nrb)
C
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nramin            : Minimum row number of A 
C integer nramax            : Maximum row number of A 
C integer nca               : Number of columns of A
C integer namax             : First dimension of A matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of A)
C integer numa(nra)         : Control vector of A matrix
C                            (number of nonzero elements of each row of A)
C integer lista(namax,nra)  : Control vector of A matrix
C                           (list of nonzero elements of each row of A)
C integer nrb               : Number of rows of B
C integer ncb               : Number of columns of B
C integer nbmax             : First dimension of B matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of B)
C integer numb(nrb)         : Control vector of B matrix
C                             (number of nonzero elements of each row of B)
C integer listbptr(nrb)     : Control vector of B matrix
C                             (pointer to the start of each row of B)
C integer listb(nbmax)      : Control vector of B matrix
C                             (list of nonzero elements of each row of B)
C integer ncmax             : First dimension of C matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of C)
C **************************** OUTPUT *********************************
C integer numc(nra)          : Control vector of C matrix
C                            (number of nonzero elements of each row of C)
C integer listc(ncmax,nra)   : Control vector of C matrix
C                            (list of nonzero elements of each row of C)
C **************************** INTERNAL *******************************
C integer ind(ncb)          : Auxiliary array to build C in sparse form
C integer nindv(ncmax)      : Auxiliary array to store indexes of nonzero
C                             matrix elements of each row of C
C *********************************************************************
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .  nca,ncb,namax,nbmax,ncmax,nramin,nramax,nrb,Node,Nodes,
     .  lista(namax,*),listb(nbmax),listc(ncmax,*),
     .  numa(*),numb(nrb),numc(*),listbptr(nrb)

C Internal variables..................................................
      integer
     .  i,il,in,j,k,kl,kn,lbmax,n,nbpmax,nraloc,nrbl

      integer, dimension(:,:), allocatable ::
     .  ind, nindv
#ifdef MPI
      integer, dimension(:), allocatable ::
     .  listbtmp
      integer
     .  MPIerror, BNode, nl, numbtmp
#endif

C Check dimensions ....................................................
      call chkdim('axb_build','nca',nca,nrb,0)

C Initialize internal variables .......................................
      nraloc = nramax - nramin + 1
      allocate(ind(ncb,nraloc))
      call memory('A','I',ncb*nraloc,'axb_build')
      allocate(nindv(ncmax,nraloc))
      call memory('A','I',ncmax*nraloc,'axb_build')
      numc(1:nraloc) = 0
      ind(1:ncb,1:nraloc) = 0
      nindv(1:ncmax,1:nraloc)=0

C Find out C control vectors ..........................................
#ifdef MPI
      call MPI_AllReduce(nbmax,nbpmax,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)
      allocate(listbtmp(nbmax))
      call memory('A','I',nbmax,'axb_build1')
#else
      nbpmax = listbptr(nrb)+numb(nrb)
#endif

C Loop over node whose control vector for b is to be distributed
      do n = 1,nrb
#ifdef MPI
        call WhichNodeOrb(n,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(n,Node,Nodes,nl)
          numbtmp = numb(nl)
          listbtmp(1:numbtmp) = listb(listbptr(nl)+1:listbptr(nl)+
     .      numbtmp)
        endif
        call MPI_Bcast(numbtmp,1,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(listbtmp,numbtmp,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
C Loop over local bands
        do i = 1,nraloc
          do in = 1,numa(i)
            k = lista(in,i)
            if (k .eq. n) then
#ifdef MPI
              do kn = 1,numbtmp
                j = listbtmp(kn)
#else
              do kn = 1,numb(n)
                j = listb(kn)
#endif
                if (ind(j,i) .eq. 0) then
                  ind(j,i) = 1
                  numc(i) = numc(i) + 1
                  nindv(numc(i),i) = j
                endif
              enddo
            endif
          enddo
        enddo
      enddo
      do i = 1,nraloc
        call chkdim ('axb_build','ncmax',ncmax,numc(i),1)
        do in = 1,numc(i)
          j = nindv(in,i)
          listc(in,i) = j
        enddo
      enddo
C ...........................
C Free local memory
      call memory('D','I',size(ind),'axb_build')
      deallocate(ind)
      call memory('D','I',size(nindv),'axb_build')
      deallocate(nindv)
#ifdef MPI
      call memory('D','I',size(listbtmp),'axb_build1')
      deallocate(listbtmp)
#endif

      return
      end

      subroutine axb_build2(nramin,nramax,nca,namax,numa,lista,
     .               nrb,ncb,nbmax,numb,listb,ncmax,numc,listc)
C ********************************************************************
C Constructs control indexes of a C matrix in sparse form,
C C being the product of A and B (also in sparse form)
C
C              C = A x B
C
C In full form: A is rectangular, and has dimension:  nra x nca
C               B is rectangular, and has dimension:  nrb x ncb
C and, as a result:
C               C is rectangular, and has dimension:  nra x ncb
C (Of course, nca must be equal to nrb)
C
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nramin            : Minimum row number of A 
C integer nramax            : Maximum row number of A 
C integer nca               : Number of columns of A
C integer namax             : First dimension of A matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of A)
C integer numa(nra)         : Control vector of A matrix
C                            (number of nonzero elements of each row of A)
C integer lista(namax,nra)  : Control vector of A matrix
C                           (list of nonzero elements of each row of A)
C integer nrb               : Number of rows of B
C integer ncb               : Number of columns of B
C integer nbmax             : First dimension of B matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of B)
C integer numb(nrb)          : Control vector of B matrix
C                            (number of nonzero elements of each row of B)
C integer listb(nbmax,nrb)   : Control vector of B matrix
C                            (list of nonzero elements of each row of B)
C integer ncmax             : First dimension of C matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of C)
C **************************** OUTPUT *********************************
C integer numc(nra)          : Control vector of C matrix
C                            (number of nonzero elements of each row of C)
C integer listc(ncmax,nra)   : Control vector of C matrix
C                            (list of nonzero elements of each row of C)
C **************************** INTERNAL *******************************
C integer ind(ncb)          : Auxiliary array to build C in sparse form
C integer nindv(ncmax)      : Auxiliary array to store indexes of nonzero
C                             matrix elements of each row of C
C *********************************************************************
      implicit none

      integer
     .  nca,ncb,namax,nbmax,ncmax,nramin,nramax,nrb,
     .  lista(namax,*),listb(nbmax,nrb),listc(ncmax,*),
     .  numa(*),numb(nrb),numc(*)

C Internal variables..................................................
      integer
     .  i,il,in,j,k,kn,nind

      integer, dimension(:), allocatable ::
     .  ind, nindv
C............................

C Check dimensions ....................................................
      call chkdim('axb_build','nca',nca,nrb,0)
C ...........................
C Initialize internal variables .......................................
      allocate(ind(ncb))
      call memory('A','I',ncb,'axb_build')
      allocate(nindv(ncmax))
      call memory('A','I',ncmax,'axb_build')
      nind=0
      do i = 1,ncb
        ind(i) = 0
      enddo
      do i = 1,ncmax
        nindv(i)=0
      enddo
C ...........................
C Find out C control vectors ..........................................
      do i = nramin,nramax
        il = i - nramin + 1
        do in = 1,numa(il)
          k = lista(in,il)
          do kn = 1,numb(k)
            j = listb(kn,k)
            if (ind(j) .eq. 0) then
              ind(j) = 1
              nind = nind+1
              nindv(nind) = j
            endif
          enddo
        enddo
        numc(il) = nind
        call chkdim ('axb_build','ncmax',ncmax,nind,1)
        do in = 1,nind
          j = nindv(in)
          nindv(in) = 0
          ind(j) = 0
          listc(in,il) = j
        enddo
        nind = 0
      enddo
C ...........................
C Free local memory
      call memory('D','I',size(ind),'axb_build')
      deallocate(ind)
      call memory('D','I',size(nindv),'axb_build')
      deallocate(nindv)

      return
      end
