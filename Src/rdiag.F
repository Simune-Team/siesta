      subroutine rdiag(H,S,n,nm,nml,maxuo,w,Z,neigvec,ierror)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C real symmetric matrices, by calling the  LAPACK routine DSYGV
C Written by G.Fabricius and J.Soler, March 1998
C Rewritten by Julian Gale, June 2004
C ************************** INPUT ******************************************
C real*8 H(nml,nm)                 : Symmetric H matrix
C real*8 S(nml,nm)                 : Symmetric S matrix
C integer n                        : Order of the generalized  system
C integer nm                       : Right hand dimension of H and S matrices
C integer nml                      : Left hand dimension of H and S matrices
C                                    which is greater than or equal to nm
C integer neigvec                  : No. of eigenvectors to calculate
C ************************** OUTPUT *****************************************
C real*8 w(nml)                    : Eigenvalues
C real*8 Z(nml,nm)                 : Eigenvectors
C integer ierror                   : Flag indicating success code for routine
C                                  :  0 = success
C                                  : -1 = repeat call as memory is increased
C                                  :  1 = fatal error
C ************************* PARALLEL ****************************************
C When running in parallel this routine now uses Scalapack to perform a
C parallel matrix diagonalisation. This requires Scalapack and Blacs to
C be installed first. Note that a one-dimensional processor grid is used
C as this simplifies the mapping on to the distributed data structure.
C Note here we call the expert driver routine form of Lapack which allows
C us the possibility of only calculating a subset of the eigenvalues. This
C is not used, but could be useful in future. Also abstol, the accuracy
C with which the eigenvalues are found, is currently set to a default.
C In future benefits could be gained from relating this value to the
C current degree of convergence of the SCF process.
C Modified by Julian Gale, November 1998
C ***************************************************************************
C
C  Modules
C
      use precision
      use fdf
      use parallel,   only : Node, Nodes, BlockSize
      use diagmemory, only : MemoryFactor
#ifdef MPI
      use mpi_siesta, only : mpi_bcast, mpi_logical, mpi_comm_world
#endif
      use alloc
      use sys,       only : die

      implicit          none

C Passed variables
      integer                 :: ierror
      integer                 :: maxuo
      integer                 :: n
      integer                 :: neigvec
      integer                 :: nm
      integer                 :: nml
      real(dp)                :: H(nml,nm)
      real(dp)                :: S(nml,nm)
      real(dp)                :: w(nml)
      real(dp)                :: Z(nml,nm)

C Local variables
      type(allocDefaults) oldDefaults
#ifdef MPI
      integer                 :: anb
      integer                 :: desch(9)
      integer                 :: MPIerror
      integer                 :: iacol
      integer                 :: iarow
      integer                 :: iceil
      integer                 :: ictxt
      integer                 :: indxg2p
      integer                 :: lwork2
      integer                 :: mq0
      integer                 :: mvfound
      integer                 :: mfound
      integer                 :: mycol
      integer                 :: myrow
      integer                 :: nn
      integer                 :: np
      integer                 :: np0
      integer                 :: npcol
      integer                 :: nprow
      integer                 :: nps
      integer                 :: nq
      integer                 :: nsygst_lwopt
      integer                 :: nsytrd_lwopt
      integer                 :: nt
      integer                 :: numroc
      integer                 :: pjlaenv
      integer                 :: sqnpc
      integer,  pointer, save :: iclustr(:)
      real(dp), pointer, save :: gap(:)
      logical                 :: BlacsOK
#endif
      character               :: jobz
      character               :: range
      integer                 :: i1
      integer                 :: i2
      integer                 :: ilaenv
      integer                 :: info
      integer                 :: liwork
      integer                 :: lwork
      integer                 :: neigok
      integer                 :: nb
      integer                 :: nz
      integer,  pointer, save :: ifail(:)
      integer,  pointer, save :: iwork(:)
      logical,           save :: DivideConquer
      logical,           save :: DCread = .true.
      logical,           save :: nullified = .false.
      logical                 :: Serial
      real(dp)                :: abstol
      real(dp)                :: one
      real(dp)                :: orfac
      real(dp)                :: scale
      real(dp)                :: vl
      real(dp)                :: vu
      real(dp), pointer, save :: work(:)

C Start time count
      call timer('rdiag',1)

C Nullify pointers
      if (.not.nullified) then
        nullify( ifail, iwork, work )
#ifdef MPI
        nullify( gap, iclustr )
#endif
        nullified = .true.
      endif

C Get old allocation defaults and set new ones
      call alloc_default( old=oldDefaults, copy=.false., shrink=.true.,
     .                    imin=1, routine='rdiag' )

C*******************************************************************************
C Setup                                                                        *
C*******************************************************************************

C Initialise error flag
      ierror = 0

C Set algorithm logicals
      Serial = (nm.eq.nml.or.Nodes.eq.1)
      if (DCread) then
        DCread = .false.
        if (Node.eq.0) then
          DivideConquer = fdf_boolean( 'DivideAndConquer', .false. )
        endif
#ifdef MPI
        call MPI_Bcast(DivideConquer,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
#endif
      endif

#ifdef MPI
C Get Blacs context and initialise Blacs grid 
      if (.not.Serial) then
        nprow = 1
        npcol = Nodes
        call blacs_get( -1, 0, ictxt )
        call blacs_gridinit( ictxt, 'C', nprow, npcol )

C Set up blacs descriptors for parallel case
        BlacsOK = .true.                                                       
        call descinit( desch, nml, nml, BlockSize, BlockSize, 0, 0, 
     .    ictxt, nml, info )
        if (info.ne.0) BlacsOK = .false.
        if (.not.BlacsOK) then
          call die('ERROR : Blacs setup has failed in rdiag!')
        endif
      endif
#endif

C Set general Lapack parameters
      if (neigvec.gt.0) then
        jobz   = 'V'
        if (neigvec.eq.nml) then
          range  = 'A'
        else
          range  = 'I'
        endif
      else
        jobz   = 'N'
        range  = 'A'
      endif
      abstol = 1.0d-8
      orfac  = 1.0d-3
      one = 1.0d0

C Calculate memory requirements
      if (Serial) then
        if (DivideConquer) then
          if (neigvec.gt.0) then
            lwork = 1 + 6*n + n*n
            liwork = 3 + 5*n
          else
            lwork = 1 + 2*n
            liwork = 1
          endif
        else
          nb = ilaenv(1,'DSYTRD','U',n,-1,-1,-1)
          lwork = max(8*n,(nb+3)*n)
          liwork = 5*n
        endif
#ifdef MPI
      else
        call blacs_gridinfo(desch(2),nprow,npcol,myrow,mycol)
        if (DivideConquer) then
          np0 = numroc(n,BlockSize,0,0,nprow)
          mq0 = numroc(n,BlockSize,0,0,npcol)
          lwork = max(BlockSize*(np0+1),3*BlockSize)
          lwork2 = BlockSize*(2*np0 + mq0 + BlockSize)
          lwork = max(lwork,lwork2)
          iarow = indxg2p(1,BlockSize,myrow,desch(7),nprow)
          iacol = indxg2p(1,BlockSize,mycol,desch(8),npcol)
          np = numroc(n,BlockSize,myrow,iarow,nprow)
          nq = numroc(n,BlockSize,mycol,iacol,npcol)
          nt = 3*n + max(BlockSize*(np+1),3*BlockSize)
          lwork2 = max(1+6*n+2*np*nq,nt) + 2*n
          lwork = max(lwork,lwork2)
          liwork = 7*n + 8*Nodes + 2
        else
          nn = max(n,BlockSize)
          np0 = numroc(nn,BlockSize,0,0,1)
          if (neigvec.gt.0) then
            mq0 = numroc(max(n,BlockSize),BlockSize,0,0,Nodes)
            lwork = 5*n + max(5*nn,np0*mq0+2*BlockSize*BlockSize) +
     .              iceil(neigvec,Nodes)*nn
            anb = pjlaenv(ictxt,3,'PDSYTTRD','L',0,0,0,0)
            sqnpc = int(sqrt(dble(Nodes)))
            nps = max(numroc(n,1,0,0,sqnpc),2*anb)
            nsytrd_lwopt = n + 2*(anb + 1)*(4*nps + 2) + 
     .                     (nps + 3)*nps
            nsygst_lwopt = (2*np + nq)*BlockSize + 
     .                     BlockSize*BlockSize
            lwork2 = max(5*n + nsytrd_lwopt,nsygst_lwopt)
            lwork = max(lwork,lwork2)
          else
            lwork = 5*n + max(5*nn,BlockSize*(np+1))
          endif
          liwork = 6*max(n,Nodes+1,4)
        endif
#endif
      endif

C Scale memory by memory factor
      lwork = nint(MemoryFactor*dble(lwork))

C Allocate workspace arrays
      call re_alloc( work,    1,lwork,   name='work'  )
      call re_alloc( iwork,   1,liwork,  name='iwork' )
      call re_alloc( ifail,   1,nml,     name='ifail'   )
#ifdef MPI
      call re_alloc( gap,     1,Nodes,   name='gap'     )
      call re_alloc( iclustr, 1,2*Nodes, name='iclustr' )
#endif

C*******************************************************************************
C Factorise overlap matrix                                                     *
C*******************************************************************************
      call timer('rdiag1',1)
      if (Serial) then
        call dpotrf('U',n,S,nml,info)
#ifdef MPI
      else
        call pdpotrf('U',n,S,1,1,desch,info)
#endif
      endif
      if (info.ne.0) then
        call die('Error in Cholesky factorisation in rdiag')
      endif
      call timer('rdiag1',2)

C*******************************************************************************
C Transform problem to standard eigenvalue problem                             *
C*******************************************************************************
      call timer('rdiag2',1)
      if (Serial) then
        call dsygst(1,'U',n,H,nml,S,nml,info)
#ifdef MPI
      else
        call pdsyngst(1,'U',n,H,1,1,desch,S,1,1,
     .                desch,scale,work,lwork,info)
#endif
      endif
      if (info.ne.0) then
        call die('Error in forward transformation in rdiag')
      endif
      call timer('rdiag2',2)

C*******************************************************************************
C Solve standard eigenvalue problem                                            *
C*******************************************************************************
      call timer('rdiag3',1)
      if (Serial) then
        if (DivideConquer) then
          call dsyevds(jobz,'U',n,H,nml,w,Z,nml,work,lwork,iwork,
     .      liwork,info)
          neigok = n
        else
          call dsyevx(jobz,range,'U',n,H,nml,vl,vu,1,neigvec,abstol,
     .                neigok,w,Z,nml,work,lwork,iwork,ifail,info)
        endif
#ifdef MPI
      else
        if (DivideConquer) then
          call pdsyevd(jobz,'U',n,H,1,1,desch,w,Z,1,1,desch,work,
     .      lwork,iwork,liwork,info)
          neigok = n
        else
          call pdsyevx(jobz,range,'U',n,H,1,1,desch,vl,vu,1,neigvec,
     .                 abstol,neigok,nz,w,orfac,Z,1,1,desch,work,
     .                 lwork,iwork,liwork,ifail,iclustr,gap,info)
        endif
#endif
      endif

C Check error flag
      if (info.ne.0) then
        ierror = 1
        if (info.lt.0) then
          call die('Illegal argument to standard eigensolver')
        elseif (info.gt.0)   then
          call die('Failure to converge standard eigenproblem')
        endif
      endif
      if (neigok.lt.neigvec) then
        call die('Insufficient eigenvalues converged in rdiag')
      endif
      call timer('rdiag3',2)

C*******************************************************************************
C Back transformation                                                          *
C*******************************************************************************
      call timer('rdiag4',1)
      if (Serial) then
        call dtrsm('Left','U','N','Non-unit',n,neigvec,one,
     .             S,nml,Z,nml)
#ifdef MPI
      else
        call pdtrsm('Left','U','N','Non-unit',n,neigvec,one,
     .              S,1,1,desch,Z,1,1,desch)
        if (scale.ne.one) then
          call dscal(n,scale,w,1)
        endif
#endif
      endif
      if (info.ne.0) then
        call die('Error in back transformation in rdiag')
      endif
      call timer('rdiag4',2)

C*******************************************************************************
C Clean up                                                                     *
C*******************************************************************************

C Deallocate workspace arrays
#ifdef MPI
      call de_alloc( iclustr, name='iclustr')
      call de_alloc( gap,     name='gap'    )
#endif
      call de_alloc( ifail,   name='ifail'  )
      call de_alloc( iwork,   name='iwork'  )
      call de_alloc( work,    name='work'   )

#ifdef MPI
C Exit Blacs
      if (.not.Serial) then
        call blacs_gridexit( ictxt )
      endif
#endif

C  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

C Stop time count
      call timer('rdiag',2)

      end
