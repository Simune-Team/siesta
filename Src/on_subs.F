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

      call timer('ctrans1',1)

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
      call timer('ctrans1',2)

      return
      end
      subroutine ctrans2s(nr,nc,nmax,ntmax,num,list,numt,Node,Nodes)
C ********************************************************************
C Finds the C maximum number of non-zero column elements, ntmax.
C Adapted from ctrans.
C Written by P.Ordejon. October'96
C Adapted by J.D.Gale. December'01
C ***************************** INPUT *********************************
C integer nr             : Number of rows of C 
C integer nc             : Number of columns of full C
C integer nmax           : First dimension of list and C, and maximum
C                           number of nonzero elements of each row of C
C integer num(nr)        : Control vector of C matrix
C                           (number of nonzero elements of each row of C)
C integer list(nmax,nr)  : Control vector of C matrix
C                          (list of nonzero elements of each row of C)
C integer Node           : Local node number
C integer Nodes          : Total number of global nodes
C **************************** OUTPUT *********************************
C integer numt(nc)       : Control vector of C matrix
C                         (number of nonzero elements of each column of C)
C integer ntmax          : maximum number of nonzero elements of each
C                           column of C
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
     .  list(nmax,nr),num(nr),numt(nc)

#ifdef MPI
      integer
     .  MPIerror, nbandspernode, nbandsmin, nbandsmax, nremainder,
     .  in,nbandsloc,nmaxl,nmaxg,nrmax,nrmin,nrloc,mul
      integer, dimension(:,:), allocatable, save :: listlocal
      integer, dimension(:)  , allocatable, save :: numlocal
#endif

C Internal variables ..................................................
      integer
     .  i,imu,mu,n,il
C ..........................
      call timer('ctrans2s',1)

C Initialize numt list ...............................................
      do i = 1,nc
        numt(i) = 0
      enddo

C Construct information for transpose of C ............................
C Need to use information from other Nodes

#ifdef MPI
C Work out local distribution of bands
      nbandspernode = ((nr-1)/Nodes)
      nremainder = nr - nbandspernode*Nodes
      nbandsmin = Node*nbandspernode + min(nremainder,Node) + 1
      nbandsmax = nbandsmin + nbandspernode - 1
      if (Node.lt.nremainder) nbandsmax = nbandsmax + 1
      nbandsloc = nbandsmax - nbandsmin + 1

C Get global maximum number of non-zero elements in a given row
      nmaxl = 0
      do i = 1,nbandsloc
        nmaxl = max(nmaxl,num(i))
      enddo
      call MPI_AllReduce(nmaxl,nmaxg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Allocate local memory
      allocate(numlocal(nbandspernode+1))
      call memory('A','I',nbandspernode+1,'ctrans2')
      allocate(listlocal(nmaxg,nbandspernode+1))
      call memory('A','I',nmaxg*(nbandspernode+1),'ctrans2')

C Loop over Nodes broadcasting num/list
      do in = 0,Nodes-1

C Work out which bands are on Node = in
        nrmin = in*nbandspernode + min(nremainder,in) + 1
        nrmax = nrmin + nbandspernode - 1
        if (in.lt.nremainder) nrmax = nrmax + 1
        nrloc = nrmax - nrmin + 1

C Copy data on broadcast node to workspace arrays
        if (Node.eq.in) then
          do i = 1,nrloc
            numlocal(i) = num(i)
            do il = 1,nmaxg
              listlocal(il,i) = list(il,i)
            enddo
          enddo
        endif

C Transfer data
        call MPI_Bcast(numlocal,nrloc,MPI_integer,in,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(listlocal(1,1),nmaxg*nrloc,MPI_integer,in,
     .    MPI_Comm_World,MPIerror)

C Use information from Node = in
        do mu = nrmin,nrmax
          mul = mu - nrmin + 1
          do imu = 1,numlocal(mul)
            i = listlocal(imu,mul)
            numt(i) = numt(i)+1
          enddo
        enddo
      enddo

C Free local memory
      call memory('D','I',size(listlocal),'ctrans2')
      deallocate(listlocal)
      call memory('D','I',size(numlocal),'ctrans2')
      deallocate(numlocal)
#else
      do mu = 1,nr
        do imu = 1,num(mu)
          i = list(imu,mu)
          numt(i) = numt(i)+1
        enddo
      enddo
#endif

C Find ntmax
      ntmax = 0
      do i = 1,nc
        ntmax = max(ntmax,numt(i))
      enddo

C ..........................
      call timer('ctrans2s',2)

      return
      end

      subroutine ctrans2(nr,nc,nmax,ntmax,
     .              num,list,numt,listt,fttoc,Node,Nodes)
C ********************************************************************
C Finds the C transpose matrix control vectors numt and listt,
C and the index vector fttoc.
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
C integer fttoc(ntmax,nc) : Map from C transpose to C indexing
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
     .  fttoc(ntmax,*),list(nmax,nr),listt(ntmax,*),
     .  num(nr),numt(nc)

#ifdef MPI
      integer
     .  MPIerror, nbandspernode, nbandsmin, nbandsmax, nremainder,
     .  in,nbandsloc,nmaxl,nmaxg,nrmax,nrmin,nrloc,mul
      integer, dimension(:,:), allocatable, save :: listlocal
      integer, dimension(:)  , allocatable, save :: numlocal
      integer, dimension(:)  , allocatable, save :: iptr
#endif

C Internal variables ..................................................
      integer
     .  i,imu,mu,n,il
C ..........................
      call timer('ctrans2',1)

C  Initialize numt list ...............................................
      do i = 1,nc
        numt(i) = 0
      enddo

C ..........................
C  Construct information for transpose of C ............................
C  Need to use information from other Nodes

#ifdef MPI
C Work out local distribution of bands
      nbandspernode = ((nr-1)/Nodes)
      nremainder = nr - nbandspernode*Nodes
      nbandsmin = Node*nbandspernode + min(nremainder,Node) + 1
      nbandsmax = nbandsmin + nbandspernode - 1
      if (Node.lt.nremainder) nbandsmax = nbandsmax + 1
      nbandsloc = nbandsmax - nbandsmin + 1

C Get global maximum number of non-zero elements in a given row
      nmaxl = 0
      do i = 1,nbandsloc
        nmaxl = max(nmaxl,num(i))
      enddo
      call MPI_AllReduce(nmaxl,nmaxg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Allocate local memory
      allocate(iptr(nc))
      call memory('A','I',nc,'ctrans2')
      allocate(numlocal(nbandspernode+1))
      call memory('A','I',nbandspernode+1,'ctrans2')
      allocate(listlocal(nmaxg,nbandspernode+1))
      call memory('A','I',nmaxg*(nbandspernode+1),'ctrans2')

C Build pointer for current Node
      do i = 1,nc
        call GlobalToLocalOrb(i,Node,Nodes,iptr(i))
      enddo

C Loop over Nodes broadcasting num/list
      do in = 0,Nodes-1

C Work out which bands are on Node = in
        nrmin = in*nbandspernode + min(nremainder,in) + 1
        nrmax = nrmin + nbandspernode - 1
        if (in.lt.nremainder) nrmax = nrmax + 1
        nrloc = nrmax - nrmin + 1

C Copy data on broadcast node to workspace arrays
        if (Node.eq.in) then
          do i = 1,nrloc
            numlocal(i) = num(i)
            do il = 1,nmaxg
              listlocal(il,i) = list(il,i)
            enddo
          enddo
        endif

C Transfer data
        call MPI_Bcast(numlocal,nrloc,MPI_integer,in,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(listlocal(1,1),nmaxg*nrloc,MPI_integer,in,
     .    MPI_Comm_World,MPIerror)

C Use information from Node = in
        do mu = nrmin,nrmax
          mul = mu - nrmin + 1
          do imu = 1,numlocal(mul)
            i = listlocal(imu,mul)
            numt(i) = numt(i)+1
            n = numt(i)
            il = iptr(i)
            if (il.gt.0) then
              listt(n,il) = mu
              fttoc(n,il) = imu
            endif
          enddo
        enddo
      enddo

C Free local memory
      call memory('D','I',size(listlocal),'ctrans2')
      deallocate(listlocal)
      call memory('D','I',size(numlocal),'ctrans2')
      deallocate(numlocal)
      call memory('D','I',size(iptr),'ctrans2')
      deallocate(iptr)
#else
      do mu = 1,nr
        do imu = 1,num(mu)
          i = list(imu,mu)
          numt(i) = numt(i)+1
          n = numt(i)
          listt(n,i) = mu
          fttoc(n,i) = imu
        enddo
      enddo
#endif
C ..........................
      call timer('ctrans2',2)

      return
      end

      subroutine axb_build1s(nramin,nramax,nca,namax,numa,lista,
     .               nrb,ncb,nbmax,numb,listbptr,listb,ncmax,
     .               numc,node,Nodes)
C ********************************************************************
C Determines the maximum number of non-zero elements in a row of a
C sparse matrix C :
C
C              C = A x B
C
C In full form: A is rectangular, and has dimension:  nra x nca
C               B is rectangular, and has dimension:  nrb x ncb
C and, as a result:
C               C is rectangular, and has dimension:  nra x ncb
C (Of course, nca must be equal to nrb)
C
C Adapted from axb_build.
C Written by P.Ordejon. October'96
C Adapted by J.D.Gale. December'01
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
C **************************** OUTPUT *********************************
C integer numc(nra)         : Control vector of C matrix
C                            (number of nonzero elements of each row of C)
C integer ncmax             : Maximum number of non-zero row elements of C
C **************************** INTERNAL *******************************
C integer ind(ncb)          : Auxiliary array to build C in sparse form
C integer nindv(ncmax,nraloc): Auxiliary array to store indexes of nonzero
C                             matrix elements of each row of C
C *********************************************************************
      use parallel
#ifdef MPI
      use mpi
#endif
      use alloc

      implicit none

      integer
     .  nca,ncb,namax,nbmax,ncmax,nramin,nramax,nrb,Node,Nodes,
     .  lista(namax,*),listb(nbmax),numa(*),numb(nrb),numc(*),
     .  listbptr(nrb)

C Internal variables..................................................
      integer
     .  i,in,j,k,kn,n,nn,nraloc,nrbloc,nrbmin,nrbmax

      integer, dimension(:), allocatable ::
     .  indl

      integer, pointer, save :: nindv(:,:)

#ifdef MPI
      integer, dimension(:), allocatable ::
     .  numbtmp, iptr
      integer, dimension(:,:), allocatable ::
     .  listbtmp
      integer
     .  MPIerror, BNode, kloc, nl, nbpl, nbpmax
#endif

      call timer('axbbuild1s',1)

C Check dimensions ....................................................
      call chkdim('axb_build1s','nca',nca,nrb,0)

C Initialize internal variables .......................................
      nraloc = nramax - nramin + 1
      allocate(indl(ncb))
      call memory('A','I',ncb,'axb_build1s')
      allocate(nindv(ncmax,nraloc))
      call memory('A','I',ncmax*nraloc,'axb_build1s')
      numc(1:nraloc) = 0
      indl(1:ncb) = 0

C Allocate initial size of nindv
      nullify(nindv)
      ncmax = 100
      call realloc(nindv,1,ncmax,1,nraloc,name='nindv')
      nindv(1:ncmax,1:nraloc) = 0

C Find out C control vectors ..........................................
#ifdef MPI
      call GetNodeOrbs(nrb,Node,Nodes,nrbloc)
      nbpl = 0
      do i = 1,nrbloc
        nbpl = max(nbpl,numb(i))
      enddo
      call MPI_AllReduce(nbpl,nbpmax,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)
      call GetNodeOrbs(nrb,0,Nodes,nrbmax)
      allocate(iptr(nrb))
      call memory('A','I',nrb,'axb_build1s')
      allocate(numbtmp(nrbmax))
      call memory('A','I',nrbmax,'axb_build1s')
      allocate(listbtmp(nbpmax,nrbmax))
      call memory('A','I',nbpmax*nrbmax,'axb_build1s')
#endif

C Loop over node whose control vector for b is to be distributed
      do nn = 1,Nodes
#ifdef MPI
        call GetNodeOrbs(nrb,nn-1,Nodes,nrbloc)
        call LocalToGlobalOrb(1,nn-1,Nodes,nrbmin)
        call LocalToGlobalOrb(nrbloc,nn-1,Nodes,nrbmax)
        if (Node.eq.nn-1) then
          numbtmp(1:nrbloc) = numb(1:nrbloc)
          do i = 1,nrbloc
            listbtmp(1:numbtmp(i),i) = listb(listbptr(i)+1:listbptr(i)+
     .      numb(i))
          enddo
        endif
        call MPI_Bcast(numbtmp,nrbloc,MPI_integer,nn-1,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(listbtmp(1,1),nbpmax*nrbloc,MPI_integer,nn-1,
     .    MPI_Comm_World,MPIerror)
        do i = 1,nrb
          call GlobalToLocalOrb(i,nn-1,Nodes,iptr(i))
        enddo
#else
        nrbmin = 1
        nrbmax = nrb
#endif
C Loop over local bands
        do i = 1,nraloc

C Expand index
          do j = 1,numc(i)
            indl(nindv(j,i)) = 1
          enddo

          do in = 1,numa(i)
            k = lista(in,i)
#ifdef MPI
            kloc = iptr(k)
            if (kloc .gt. 0) then
              do kn = 1,numbtmp(kloc)
                j = listbtmp(kn,kloc)
#else
              do kn = 1,numb(k)
                j = listb(listbptr(k)+kn)
#endif
                if (indl(j) .eq. 0) then
                  indl(j) = 1
                  numc(i) = numc(i) + 1
                  if (numc(i).gt.ncmax) then
                    ncmax = ncmax + 100
                    call realloc(nindv,1,ncmax,1,nraloc,name='nindv')
                  endif
                  nindv(numc(i),i) = j
                endif
              enddo
#ifdef MPI
            endif
#endif
          enddo

C Compress index and re-zero indl
          do j = 1,numc(i)
            indl(nindv(j,i)) = 0
          enddo

        enddo
      enddo
      ncmax = 0
      do i = 1,nraloc
        ncmax = max(ncmax,numc(i))
      enddo
C ...........................

C Free local memory
      call dealloc( nindv, name='nindv' )
      call memory('D','I',size(indl),'axb_build1s')
      deallocate(indl)
#ifdef MPI
      call memory('D','I',size(iptr),'axb_build1s')
      deallocate(iptr)
      call memory('D','I',size(numbtmp),'axb_build1s')
      deallocate(numbtmp)
      call memory('D','I',size(listbtmp),'axb_build1s')
      deallocate(listbtmp)
#endif

      call timer('axbbuild1s',2)

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
     .  lista(namax,*),listb(nbmax),listc(ncmax,*),numa(*),
     .  numb(nrb),numc(*),listbptr(nrb)

C Internal variables..................................................
      integer
     .  i,in,j,k,kn,n,nn,nraloc,nrbloc,nrbmin,nrbmax

      integer, dimension(:), allocatable ::
     .  indl

      integer, dimension(:,:), allocatable ::
     .  nindv
#ifdef MPI
      integer, dimension(:), allocatable ::
     .  numbtmp, iptr
      integer, dimension(:,:), allocatable ::
     .  listbtmp
      integer
     .  MPIerror, BNode, kloc, nl, nbpl, nbpmax
#endif

      call timer('axb_build1',1)

C Check dimensions ....................................................
      call chkdim('axb_build1','nca',nca,nrb,0)

C Initialize internal variables .......................................
      nraloc = nramax - nramin + 1
      allocate(indl(ncb))
      call memory('A','I',ncb,'axb_build')
      allocate(nindv(ncmax,nraloc))
      call memory('A','I',ncmax*nraloc,'axb_build')
      numc(1:nraloc) = 0
      indl(1:ncb) = 0
      nindv(1:ncmax,1:nraloc)=0

C Find out C control vectors ..........................................
#ifdef MPI
      call GetNodeOrbs(nrb,Node,Nodes,nrbloc)
      nbpl = 0
      do i = 1,nrbloc
        nbpl = max(nbpl,numb(i))
      enddo
      call MPI_AllReduce(nbpl,nbpmax,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)
      call GetNodeOrbs(nrb,0,Nodes,nrbmax)
      allocate(iptr(nrb))
      call memory('A','I',nrb,'axb_build1')
      allocate(numbtmp(nrbmax))
      call memory('A','I',nrbmax,'axb_build1')
      allocate(listbtmp(nbpmax,nrbmax))
      call memory('A','I',nbpmax*nrbmax,'axb_build1')
#endif

C Loop over node whose control vector for b is to be distributed
      do nn = 1,Nodes
#ifdef MPI
        call GetNodeOrbs(nrb,nn-1,Nodes,nrbloc)
        call LocalToGlobalOrb(1,nn-1,Nodes,nrbmin)
        call LocalToGlobalOrb(nrbloc,nn-1,Nodes,nrbmax)
        if (Node.eq.nn-1) then
          numbtmp(1:nrbloc) = numb(1:nrbloc)
          do i = 1,nrbloc
            listbtmp(1:numbtmp(i),i) = listb(listbptr(i)+1:listbptr(i)+
     .      numb(i))
          enddo
        endif
        call MPI_Bcast(numbtmp,nrbloc,MPI_integer,nn-1,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(listbtmp(1,1),nbpmax*nrbloc,MPI_integer,nn-1,
     .    MPI_Comm_World,MPIerror)
        do i = 1,nrb
          call GlobalToLocalOrb(i,nn-1,Nodes,iptr(i))
        enddo
#else
        nrbmin = 1
        nrbmax = nrb
#endif
C Loop over local bands
        do i = 1,nraloc

C Expand index
          do j = 1,numc(i)
            indl(nindv(j,i)) = 1
          enddo

          do in = 1,numa(i)
            k = lista(in,i)
#ifdef MPI
            kloc = iptr(k)
            if (kloc .gt. 0) then
              do kn = 1,numbtmp(kloc)
                j = listbtmp(kn,kloc)
#else
              do kn = 1,numb(k)
                j = listb(listbptr(k)+kn)
#endif
                if (indl(j) .eq. 0) then
                  indl(j) = 1
                  numc(i) = numc(i) + 1
                  nindv(numc(i),i) = j
                endif
              enddo
#ifdef MPI
            endif
#endif
          enddo

C Compress index and re-zero indl
          do j = 1,numc(i)
            indl(nindv(j,i)) = 0
          enddo

        enddo
      enddo
      do i = 1,nraloc
        call chkdim ('axb_build1','ncmax',ncmax,numc(i),1)
        do in = 1,numc(i)
          j = nindv(in,i)
          listc(in,i) = j
        enddo
      enddo
C ...........................

C Free local memory
      call memory('D','I',size(indl),'axb_build1')
      deallocate(indl)
      call memory('D','I',size(nindv),'axb_build1')
      deallocate(nindv)
#ifdef MPI
      call memory('D','I',size(iptr),'axb_build1')
      deallocate(iptr)
      call memory('D','I',size(numbtmp),'axb_build1')
      deallocate(numbtmp)
      call memory('D','I',size(listbtmp),'axb_build1')
      deallocate(listbtmp)
#endif

      call timer('axb_build1',2)

      return
      end

      subroutine axb_build2s(nramin,nramax,nca,namax,numa,lista,
     .               nrb,ncb,nbmax,numb,listb,ncmax,numc)
C ********************************************************************
C Finds the maximum dimension for the number of non-zero values in
C any row of the matrix C.
C
C              C = A x B
C
C In full form: A is rectangular, and has dimension:  nra x nca
C               B is rectangular, and has dimension:  nrb x ncb
C and, as a result:
C               C is rectangular, and has dimension:  nra x ncb
C (Of course, nca must be equal to nrb)
C
C Adapted from axb_build.
C Written by P.Ordejon. October'96
C Adapted by J.D.Gale. December'01
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
C integer numb(nrb)         : Control vector of B matrix
C                            (number of nonzero elements of each row of B)
C integer listb(nbmax,nrb)  : Control vector of B matrix
C                            (list of nonzero elements of each row of B)
C **************************** OUTPUT *********************************
C integer numc(nra)         : Control vector of C matrix
C                            (number of nonzero elements of each row of C)
C integer ncmax             : Maximum value of numc( )
C **************************** INTERNAL *******************************
C integer ind(ncb)          : Auxiliary array to build C in sparse form
C integer nindv(ncb)        : Auxiliary array to store indexes of nonzero
C                             matrix elements of each row of C
C *********************************************************************
      implicit none

      integer
     .  nca,ncb,namax,nbmax,ncmax,nramin,nramax,nrb,
     .  lista(namax,*),listb(nbmax,nrb),
     .  numa(*),numb(nrb),numc(*)

C Internal variables..................................................
      integer
     .  i,il,in,j,k,kn,nind

      integer, dimension(:), allocatable ::
     .  ind, nindv

      call timer('axbbuild2s',1)

C Check dimensions
      call chkdim('axb_build2s','nca',nca,nrb,0)

C Initialise return parameter
      ncmax = 0

C Initialize internal variables
      allocate(ind(ncb))
      call memory('A','I',ncb,'axb_build2s')
      allocate(nindv(ncb))
      call memory('A','I',ncb,'axb_build2s')
      nind=0
      do i = 1,ncb
        ind(i) = 0
        nindv(i) = 0
      enddo

C Find out C control vectors
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
        ncmax = max(ncmax,nind)
        do in = 1,nind
          j = nindv(in)
          nindv(in) = 0
          ind(j) = 0
        enddo
        nind = 0
      enddo

C Free local memory
      call memory('D','I',size(ind),'axb_build2s')
      deallocate(ind)
      call memory('D','I',size(nindv),'axb_build2s')
      deallocate(nindv)

      call timer('axbbuild2s',2)

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

      call timer('axb_build2',1)

C Check dimensions ....................................................
      call chkdim('axb_build2','nca',nca,nrb,0)
C ...........................
C Initialize internal variables .......................................
      allocate(ind(ncb))
      call memory('A','I',ncb,'axb_build2')
      allocate(nindv(ncmax))
      call memory('A','I',ncmax,'axb_build2')
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
        call chkdim ('axb_build2','ncmax',ncmax,nind,1)
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
      call memory('D','I',size(ind),'axb_build2')
      deallocate(ind)
      call memory('D','I',size(nindv),'axb_build2')
      deallocate(nindv)

      call timer('axb_build2',2)

      return
      end

      subroutine ind_gf(nr,nrloc,ncmin,ncmax,nc,maxnc,nfmax,
     .              numc,listc,numf,listf,indgf,Node,Nodes)
C ********************************************************************
C Maps the F matrix into C matrix
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nr             : Number of rows of C (columns of full F)
C integer nrloc          : Number of rows of C on local Node
C integer ncmin          : First column of full C (rows of F) stored 
C                        : locally.
C integer ncmax          : Last column of full C (rows of F) stored 
C                        : locally.
C integer nc             : Last column of full C globally
C integer maxnc          : First dimension of listc and C, and maximum
C                           number of nonzero elements of each row of C
C integer nfmax          : First dimension of listf and F, and maximum
C                           number of nonzero elements of each row of F
C integer numc(nr)       : Control vector of C matrix
C                           (number of nonzero elements of each row of C)
C integer listc(ncmax,nr): Control vector of C matrix
C                          (list of nonzero elements of each row of C)
C integer numf(nc)       : Control vector of F matrix
C                           (number of nonzero elements of each row of F)
C integer listf(nfmax,nc): Control vector of F matrix
C                          (list of nonzero elements of each row of F)
C **************************** OUTPUT *********************************
C integer indgf(ncmax,nrloc) : Map from F to C
C                    indgf(i,j) is the column index (in sparse notation)
C                    of the element of F corresponding to the element
C                    (i,j) of C (in sparse notation)
C *********************************************************************

      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .  maxnc,nc,ncmin,ncmax,nfmax,nn,nr,nrloc,Node,Nodes,
     .  indgf(maxnc,nrloc),listc(maxnc,nr),listf(nfmax,*),
     .  numc(nr),numf(*)

      integer
     .  i,iloc,imu,jk,mu,mug,mmu
#ifdef MPI
      integer
     .  MPIerror, nbl, nbg, nbmin, nbmax, nbandspernode,
     .  nbandsloc, nbandsmax, nbandsmin, nfmaxl, nfmaxg,
     .  nremainder

      integer, allocatable, save ::
     .  numfl(:), listfl(:,:)
#endif

      call timer('ind_gf',1)

C Initialise array
      indgf(1:maxnc,1:nrloc) = 0

#ifdef MPI
C Find maximum number of local bands
      nbandspernode = ((nc-1)/Nodes)
      nremainder = nc - nbandspernode*Nodes
      nbl = ncmax - ncmin + 1
      call MPI_AllReduce(nbl,nbg,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Work out local distribution of bands
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

C Allocate workspace array
      allocate(numfl(nbg))
      call memory('A','I',nbg,'ind_gf')
      allocate(listfl(nfmaxg,nbg))
      call memory('A','I',nfmaxg*nbg,'ind_gf')

      do nn = 1,Nodes

C Get local size on broadcast node
        nbmin = (nn-1)*nbandspernode + min(nremainder,(nn-1)) + 1
        nbmax = nbmin + nbandspernode - 1
        if ((nn-1).lt.nremainder) nbmax = nbmax + 1
        nbl = (nbmax-nbmin+1)

C Transfer local values on Node nn into broadcast arrays
        if (Node.eq.nn-1) then
          numfl(1:nbl) = numf(1:nbl)
          listfl(1:nfmaxg,1:nbl) = listf(1:nfmaxg,1:nbl)
        endif

C Broadcast indexing arrays
        call MPI_Bcast(numfl,nbl,MPI_integer,
     .    nn-1,MPI_Comm_World,MPIerror)
        call MPI_Bcast(listfl(1,1),nfmaxg*nbl,MPI_integer,
     .    nn-1,MPI_Comm_World,MPIerror)

        do mu=1,nrloc
          call LocalToGlobalOrb(mu,Node,Nodes,mug)
          do imu=1,numc(mug)
            i=listc(imu,mug)
            if (i.ge.nbmin.and.i.le.nbmax) then
              iloc = i - nbmin + 1
              do jk=1,numfl(iloc)
                mmu=listfl(jk,iloc)
                if (mmu .eq. mug) indgf(imu,mu)=jk
              enddo
            endif
          enddo
        enddo

      enddo

C Free local workspace arrays
      call memory('D','I',size(numfl),'ind_gf')
      deallocate(numfl)
      call memory('D','I',size(listfl),'ind_gf')
      deallocate(listfl)
#else
      do mu=1,nr
        do imu=1,numc(mu)
          i=listc(imu,mu)
          do jk=1,numf(i)
            mmu=listf(jk,i)
            if (mmu .eq. mu) indgf(imu,mu)=jk
          enddo
        enddo
      enddo
#endif

      call timer('ind_gf',2)

      return
      end
