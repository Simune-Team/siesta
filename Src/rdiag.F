      subroutine rdiag(h,s,n,nm,nml,maxuo,w,z,ierror)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C real symmetric matrices, by calling the  LAPACK routine DSYGV
C Writen by G.Fabricius and J.Soler, March 1998
C ************************** INPUT ******************************************
C real*8 h(nml,nm)                 : Symmetric H matrix
C real*8 s(nml,nm)                 : Symmetric S matrix
C integer n                        : Order of  the generalized  system
C integer nm                       : Right hand dimension of H and S matrices
C integer nml                      : Left hand dimension of H and S matrices
C                                    which is greater than or equal to nm
C ************************** OUTPUT *****************************************
C real*8 w(nml)                    : Eigenvalues
C real*8 z(nml,maxuo)              : Eigenvectros
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
#ifdef MPI
      use mpi_siesta
      use parallel
#endif
      use diagmemory
      use alloc

      implicit          none

      type(allocDefaults) oldDefaults
#ifdef MPI

      integer                :: desch(9)
      integer                :: descs(9)
      integer                :: descz(9)
      integer                :: MPIerror
      integer                :: nprow
      integer                :: npcol
      integer                :: ictxt
      integer                :: mvfound
      integer                :: mfound
      integer                :: nmmax
      integer, pointer, save :: iclustr(:)
      integer, pointer, save :: ifail(:)
      real*8                 :: abstol
      real*8                 :: orfac
      real*8,  pointer, save :: gap(:)
      logical                :: BlacsOK
      logical                :: Serial
#endif
      character              :: jobz
      character              :: range
      character              :: uplo
      integer                :: i1
      integer                :: i2
      integer                :: ierror
      integer                :: info
      integer                :: liwork
      integer                :: lwork
      integer                :: maxuo
      integer                :: n
      integer                :: nm
      integer                :: nml
      integer                :: Node
      integer                :: Nodes
      integer, pointer, save :: iwork(:)
      logical,          save :: DivideConquer
      logical,          save :: DCread
      logical,          save :: nullified = .false.
      real*8                 :: h(nml,nm)
      real*8                 :: s(nml,nm)
      real*8                 :: w(nml)
      real*8                 :: z(nml,maxuo)
      real*8,  pointer, save :: work(:)

      data DCread /.false./
C ......................

C Start time count
      call timer('rdiag',1)

C Nullify pointers
      if (.not.nullified) then
        nullify( iwork )
        nullify( work )
#ifdef MPI
        nullify( gap, iclustr, ifail )
#endif
        nullified = .true.
      endif

C Get old allocation defaults and set new ones
      call alloc_default( old=oldDefaults,
     .                    copy=.false., shrink=.true.,
     .                    imin=1, routine='rdiag' )

C Initialise error flag
      ierror = 0

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
      Serial = (nm.eq.nml.or.Nodes.eq.1)
#else
      Node = 0
      Nodes = 1
#endif

C Set general Lapack parameters
      jobz='V'
      uplo='U'
      range='A'

C Set divide and conquer option - serial only
      if (Nodes.eq.1) then
        if (.not. DCread) then
          DivideConquer = fdf_boolean( 'DivideAndConquer', .false. )
          DCread=.true.
        endif
      else
        DivideConquer = .false.
      endif

#ifdef MPI
      if (.not.Serial) then
C*****************************
C  Parallel diagonalisation  *
C*****************************
C Accuracy for eigenvalues
        abstol = 1.0d-8
C Accuracy for orthogonalisation
        orfac = 1.0d-3

C Calculate sizes of workspace arrays and add margin of error
        nmmax = (nml/Nodes) + BlockSize
        lwork = 5*nml+max(5*nml,nmmax*nml+2*BlockSize*BlockSize)+
     .          ((nml/Nodes)+1)*nml
        liwork = 6*max(nml,Nodes+1,4)
        lwork = nint(MemoryFactor*dble(lwork))

C Allocate workspace arrays
        call re_alloc( work,    1,lwork,   name='work'  )
        call re_alloc( iwork,   1,liwork,  name='iwork' )
        call re_alloc( gap,     1,Nodes,   name='gap'     )
        call re_alloc( iclustr, 1,2*Nodes, name='iclustr' )
        call re_alloc( ifail,   1,nml,     name='ifail'   )

C Get Blacs context and initialise Blacs grid 
        nprow = 1
        npcol = Nodes
        call blacs_get( -1, 0, ictxt )
        call blacs_gridinit( ictxt, 'C', nprow, npcol )

C Distribute matrices over grid
        BlacsOK = .true.
        call descinit( desch, nml, nml, BlockSize, BlockSize, 0, 0, 
     .    ictxt, nml, info )
        if (info.ne.0) BlacsOK = .false.
        call descinit( descs, nml, nml, BlockSize, BlockSize, 0, 0, 
     .    ictxt, nml, info )
        if (info.ne.0) BlacsOK = .false.
        call descinit( descz, nml, nml, BlockSize, BlockSize, 0, 0, 
     .    ictxt, nml, info )
        if (info.ne.0) BlacsOK = .false.
        if (.not.BlacsOK) then
        write(6,'('' ERROR : Blacs has failed to initialise in rdiag!''
     .    )')
        call MPI_Abort(MPI_Comm_World,1,MPIerror)
        endif

C Call Scalapack to diagonalise matrix
#ifdef CRAY
        call pssygvx
#else
        call pdsygvx
#endif
     .     (1,jobz,range,uplo,nml,h,1,1,desch,s,1,1,descs,
     .      0.0d0,0.0d0,0,0,abstol,mfound,mvfound,w,orfac,z,1,1,descz,
     .      work,lwork,iwork,liwork,ifail,iclustr,gap,info)

C Exit Blacs
        call blacs_gridexit( ictxt )

C Deallocate workspace arrays
        call de_alloc( iwork,   name='iwork'  )
        call de_alloc( work,    name='work'   )
        call de_alloc( ifail,   name='ifail'  )
        call de_alloc( iclustr, name='iclustr')
        call de_alloc( gap,     name='gap'    )
      else
#endif
C***************************
C  Serial diagonalisation  *
C***************************
        if (DivideConquer) then

C Set size of workspace arrays
          lwork = 1 + 6*nm + 2*nm*nm
          liwork = 3 + 5*nm

C Allocate workspace arrays
          call re_alloc( work,  1,lwork,  name='work'  )
          call re_alloc( iwork, 1,liwork, name='iwork' )

C Call lapack
#ifdef CRAY
          call ssygvd
#else
          call dsygvd
#endif
     .      (1,jobz,uplo,n,h,nm,s,nm,w,work,lwork,iwork,liwork,info)

C Deallocate workspace arrays
          call de_alloc( iwork,   name='iwork'  )
          call de_alloc( work,    name='work'   )

        else

C Set size of workspace array
          lwork = 5*nm

C Allocate workspace array
          call re_alloc( work, 1,lwork, name='work' )

C Call lapack
#ifdef CRAY
          call ssygv
#else
          call dsygv
#endif
     .      (1,jobz,uplo,n,h,nm,s,nm,w,work,lwork,info)

C Deallocate workspace array
          call de_alloc( work, name='work' )
        endif

C Copy eigenvectors into appropriate array
        do i1 = 1,n
          do i2 = 1,n
            z(i2,i1) = h(i2,i1)
          enddo
        enddo
#ifdef MPI
      endif
#endif

C Check error flag
      if (info .ne. 0) then
        ierror = 1
        if (Node.eq.0) then
          if (info.lt.0) then
            write(6,*) 'rdiag: ERROR in routine dsygv:'
            write(6,*) -info,'-th argument had an illegal value'
          elseif (info.le.n)   then
#ifdef MPI
            if (TryMemoryIncrease) then
              if (MemoryScale.eq.1.0d0) then
                ierror = -2
              else
                ierror = -1
                MemoryFactor = MemoryScale*MemoryFactor
                write(6,*) 'Diagonalisation failed to converged'
                write(6,*) 'Memory factor increased to ',MemoryFactor
              endif
            else
              write(6,*) 'DSYGV failed to converge. ',info,
     .          ' off-diagonal elements of an intermediate '
              write(6,*) 'tridiagonal form did not converge to zero'
            endif
#else
            write(6,*) 'DSYGV failed to converge. ',info,
     .        ' off-diagonal elements of an intermediate '
            write(6,*) 'tridiagonal form did not converge to zero'
#endif
          else
            write(6,*)'The leading minor of order ',info-n,
     .        ' of B is not positive definite.'
            write(6,*)'The factorization of B could not be completed',
     .        ' and no eigenvalues or eigenvectors were computed.'
          endif
        else
          if (info.le.n.and.info.gt.0.and.TryMemoryIncrease) then
            if (MemoryScale.eq.1.0d0) then
              ierror = -2
            else
              ierror = -1
              MemoryFactor = MemoryScale*MemoryFactor
            endif
          endif
        endif
      endif

C  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

C Stop time count
      call timer('rdiag',2)
      end
