      subroutine cdiag(h_2,nh,s_2,ns,n,w,z_2,nz,nevec,ierror)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C complex hermitian matrices, by calling the  LAPACK routine ZHEGV
C Writen by G.Fabricius, J.Soler, and J.Gale. Dec. 2000
C ************************** INPUT ******************************************
C complex*16 h(nh,n)   : Hermitian H matrix (destroyed on output)
C integer    nh        : First dimension of H matrix
C complex*16 s(ns,n)   : Hermitian S matrix (destroyed on output)
C integer    ns        : First dimension of S matrix
C integer    n         : Order of the generalized system
C integer    nz        : First dimension of Z matrix
C integer    nevec     : Number of lowest eigenvectors required
C                      : (if nevec.ne.0 all eigenvectors are returned)
C ************************** OUTPUT *****************************************
C real*8     w(nh)     : Eigenvalues
C complex*16 z(nz,n)   : Eigenvectors (if nevec=0, they are not calculated,
C                        but the array is still required)
C integer    ierror    : Flag indicating success code for routine
C                      :  0 = success
C                      :  1 = fatal error
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
C Error flag and check to see whether parallel form is necessary or
C not added by Julian Gale, November 1999
C ***************************************************************************
C
C Modules
      use precision, only: dp
      use fdf
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
      use parallel
#endif
      use diagmemory, only: memoryfactor, memoryscale, trymemoryincrease
      use alloc

      implicit          none

C Arguments
      integer           ierror,n,nevec,nh,ns,nz
      double precision  w(nh)
!
!     AG: To avoid type checking errors at runtime
!
      real*8 h_2(2,nh,n), s_2(2,ns,n), z_2(2,nz,n)

C Internal varables and arrays
      character               :: jobz
      character               :: uplo
      integer                 :: i1
      integer                 :: i2
      integer                 :: info
      integer                 :: lwork
      integer                 :: liwork
      integer                 :: lrwork
      integer                 :: Node
      integer                 :: Nodes
      integer,  pointer, save :: iwork(:)
      logical,           save :: nullified = .false.
      logical,           save :: DivideConquer
      logical,           save :: DCread
      real*8,   pointer, save :: rwork(:)
      real*8,   pointer, save :: work(:,:)
      type(allocDefaults) oldDefaults
#ifdef MPI
      integer                 :: desch(9)
      integer                 :: descs(9)
      integer                 :: descz(9)
      integer                 :: MPIerror
      integer                 :: nprow
      integer                 :: npcol
      integer                 :: ictxt
      integer                 :: mfound
      integer                 :: mvfound
      integer,  pointer, save :: iclustr(:)
      integer,  pointer, save :: ifail(:)
      real*8                  :: abstol
      real*8                  :: orfac
      real*8,   pointer, save :: gap(:)
      complex*16              :: dwork(1)
      real*8                  :: drwork(1)
      integer                 :: diwork(1)
      logical                 :: BlacsOK
      logical                 :: Serial
      character               :: range
#endif

      data DCread /.false./
C ......................

C Start time count
      call timer('cdiag',1)

C Nullify pointers
      if (.not.nullified) then
        nullify( iwork )
        nullify( rwork )
        nullify( work )
#ifdef MPI
        nullify( gap, iclustr, ifail )
#endif
        nullified = .true.
      endif

C Get old allocation defaults and set new ones
      call alloc_default( old=oldDefaults,
     .                    copy=.false., shrink=.true.,
     .                    imin=1, routine='cdiag' )
      
C Initialise error code
      ierror = 0

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
      Serial = (n.eq.nh.or.Nodes.eq.1)
#else
      Node = 0
      Nodes = 1
#endif

C Set general Lapack parameters
      if (nevec.le.0) then
        jobz='N'
      else
        jobz='V'
      endif
      uplo='U'

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
        range='A'
C Accuracy for eigenvalues
        abstol=1.0d-8
C Accuracy for orthogonalisation
        orfac=1.0d-3
C Allocate workspace arrays
        call re_alloc( gap,     1,Nodes,   name='gap'     )
        call re_alloc( iclustr, 1,2*Nodes, name='iclustr' )
        call re_alloc( ifail,   1,nh,      name='ifail'   )

C Get Blacs context and initialise Blacs grid
        nprow = 1
        npcol = Nodes
        call blacs_get( -1, 0, ictxt )
        call blacs_gridinit( ictxt, 'C', nprow, npcol )

C Distribute matrices over grid
        BlacsOK = .true.
        call descinit( desch, nh, nh, Blocksize, Blocksize, 0, 0,
     .    ictxt, nh, info )
        if (info.ne.0) BlacsOK = .false.
        call descinit( descs, nh, nh, Blocksize, Blocksize, 0, 0,
     .    ictxt, nh, info )
        if (info.ne.0) BlacsOK = .false.
        call descinit( descz, nh, nh, Blocksize, Blocksize, 0, 0,
     .    ictxt, nh, info )
        if (info.ne.0) BlacsOK = .false.
        if (.not.BlacsOK) then
        write(6,'('' ERROR : Blacs has failed to initialise in cdiag!''
     .    )')
        call MPI_Abort(MPI_Comm_World,1,MPIerror)
        endif

C Call Scalapack to get sizes
        lwork = -1
        lrwork = -1
        liwork = -1
#ifdef CRAY
        call pchegvx
#else
        call pzhegvx
#endif
     .   (1,jobz,range,uplo,nh,h_2,1,1,desch,s_2,1,1,descs,0.0d0,0.0d0,
     .    0,0,abstol,mfound,mvfound,w,orfac,z_2,1,1,descz,dwork,lwork,
     .    drwork,lrwork,diwork,liwork,ifail,iclustr,gap,info)
        lwork=nint(MemoryFactor*real(dwork(1)))
        lrwork=nint(MemoryFactor*drwork(1))
        liwork=diwork(1)

C Allocate arrays to correct sizes
        call re_alloc( work, 1,2, 1,lwork,  name='work'  )
        call re_alloc( rwork,     1,lrwork, name='rwork' )
        call re_alloc( iwork,     1,liwork, name='iwork' )

C Repeated here as some machines corrupt these parameters
        if (nevec.le.0) then
          jobz='N'
        else
          jobz='V'
        endif
        uplo='U'
        range='A'

C Call Scalapack to run
#ifdef CRAY
        call pchegvx
#else
        call pzhegvx
#endif
     .   (1,jobz,range,uplo,nh,h_2,1,1,desch,s_2,1,1,descs,0.0d0,0.0d0,
     .    0,0,abstol,mfound,mvfound,w,orfac,z_2,1,1,descz,work,lwork,
     .    rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)

C Exit Blacs
        call blacs_gridexit( ictxt )

C Deallocate workspace arrays
        call de_alloc( iwork,   name='iwork'  )
        call de_alloc( rwork,   name='rwork'  )
        call de_alloc( work,    name='work'   )
        call de_alloc( ifail,   name='ifail'  )
        call de_alloc( iclustr, name='iclustr')
        call de_alloc( gap,     name='gap'    )
      else
#endif
C***************************
C  Serial diagonalisation  *
C***************************
C Serial: call lapack routine chegv/zhegv
C Matrix z used for auxiliary space 

        if (DivideConquer) then

C Set size of workspace arrays
          lwork = 2*n + 2*n*n
          lrwork = 1 + 5*n + 2*n*n
          liwork = 3 + 5*n

C Allocate workspace arrays
          call re_alloc( work, 1,2, 1,lwork,  name='work'  )
          call re_alloc( rwork,     1,lrwork, name='rwork' )
          call re_alloc( iwork,     1,liwork, name='iwork' )

C Call lapack
#ifdef CRAY
          call chegvd
#else
          call zhegvd
#endif
     .      (1,jobz,uplo,n,h_2,nh,s_2,ns,w,work,lwork,rwork,lrwork,
     .       iwork,liwork,info)

C Deallocate workspace arrays
          call de_alloc( iwork,   name='iwork'  )
          call de_alloc( rwork,   name='rwork'  )
          call de_alloc( work,    name='work'   )

        else

C Set size of workspace arrays
          lwork = nh*n
          lrwork = 3*n - 2

C Allocate workspace array
          call re_alloc( rwork, 1,lrwork, name='rwork' )

#ifdef CRAY
          call chegv
#else
          call zhegv
#endif
     .      (1,jobz,uplo,n,h_2,nh,s_2,ns,w,z_2,lwork,rwork,info)

C Deallocate workspace array
          call de_alloc( rwork, name='rwork' )
        endif

C Move eigenvectors to correct storage location
        do i2 = 1,n
          do i1 = 1,nh
            z_2(1,i1,i2)=h_2(1,i1,i2)
            z_2(2,i1,i2)=h_2(2,i1,i2)
          enddo
        enddo

#ifdef MPI
      endif
#endif

      if (info .ne. 0) then
        ierror = 1
        if (Node.eq.0) then
          if (info.lt.0) then
            write(6,*) 'cdiag: ERROR in routine zhegv:'
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
              write(6,*) 'ZHEGV failed to converge. ',info,
     .          ' off-diagonal elements of an intermediate '
              write(6,*) 'tridiagonal form did not converge to zero'
            endif
#else
            write(6,*) 'ZHEGV failed to converge. ',info,
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

C  Stop time count
      call timer('cdiag',2)
      end
