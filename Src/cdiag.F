      subroutine cdiag(h,nh,s,ns,n,w,z,nz,fv,ierror)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C complex hermitian matrices, by calling the  LAPACK routine ZHEGV
C Writen by G.Fabricius and J.Soler, August 1998
C ************************** INPUT ******************************************
C complex*16 h(nh,n)   : Hermitian H matrix (destroyed on output)
C integer    nh        : First dimension of H matrix
C complex*16 s(ns,n)   : Hermitian S matrix (destroyed on output)
C integer    ns        : First dimension of S matrix
C integer    n         : Order of the generalized system
C integer    nz        : First dimension of Z matrix
C ************************** OUTPUT *****************************************
C real*8     w(nh)     : Eigenvalues
C complex*16 z(nz,n)   : Eigenvectors
C integer    ierror    : Flag indicating success code for routine
C                      :  0 = success
C                      : -1 = repeat call as memory has been increased
C                      :  1 = fatal error
C ************************* AUXILIARY ***************************************
C real*8     fv(5*nh)  : Auxiliary storage array
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
C  Modules
C
      use precision
#ifdef MPI
      use mpi
      use parallel
#endif
      use diagmemory

      implicit          none

      integer           i1,i2,info,lwork,n,nh,ns,nz,Node,Nodes,ierror
      character         jobz,uplo
      double precision  w(nh), fv(5*nh)
      double complex    h(nh,n), s(ns,n), z(nz,n)

#ifdef MPI
      integer           MPIerror, desch(9), descs(9), descz(9)
      integer           nprow, npcol, ictxt, lrwork, liwork,
     .                  mfound, mvfound
      integer,          dimension(:), allocatable, save :: 
     .                  iwork,iclustr,ifail
      double precision  abstol,orfac
      double precision, dimension(:), allocatable, save :: 
     .                  gap, rwork
      double complex,   dimension(:), allocatable, save :: 
     .                  work
      double complex    dwork(1)
      double precision  drwork(1)
      integer           diwork(1)
      logical           BlacsOK, Serial
      character         range
#endif
C ......................

C Start time count
      call timer('cdiag',1)

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
      jobz='V'
      uplo='U'
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
        allocate(iclustr(2*Nodes))
        call memory('A','I',2*Nodes,'cdiag')
        allocate(gap(Nodes))
        call memory('A','D',Nodes,'cdiag')
        allocate(ifail(nh))
        call memory('A','I',nh,'cdiag')

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
     .   (1,jobz,range,uplo,nh,h,1,1,desch,s,1,1,descs,0.0d0,0.0d0,
     .    0,0,abstol,mfound,mvfound,w,orfac,z,1,1,descz,dwork,lwork,
     .    drwork,lrwork,diwork,liwork,ifail,iclustr,gap,info)
        lwork=nint(MemoryFactor*real(dwork(1)))
        lrwork=nint(MemoryFactor*drwork(1))
        liwork=diwork(1)

C Allocate arrays to correct sizes
        allocate(work(lwork))
        call memory('A','Z',lwork,'cdiag')
        allocate(rwork(lrwork))
        call memory('A','D',lrwork,'cdiag')
        allocate(iwork(liwork))
        call memory('A','I',liwork,'cdiag')
C Repeated here as some machines corrupt these parameters
        jobz='V'
        uplo='U'
        range='A'

C Call Scalapack to run
#ifdef CRAY
        call pchegvx
#else
        call pzhegvx
#endif
     .   (1,jobz,range,uplo,nh,h,1,1,desch,s,1,1,descs,0.0d0,0.0d0,
     .    0,0,abstol,mfound,mvfound,w,orfac,z,1,1,descz,work,lwork,
     .    rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)

C Exit Blacs
        call blacs_gridexit( ictxt )

C Deallocate workspace arrays
        call memory('D','Z',size(work),'cdiag')
        deallocate(work)
        call memory('D','D',size(rwork),'cdiag')
        deallocate(rwork)
        call memory('D','I',size(iwork),'cdiag')
        deallocate(iwork)
        call memory('D','I',size(iclustr),'cdiag')
        deallocate(iclustr)
        call memory('D','D',size(gap),'cdiag')
        deallocate(gap)
        call memory('D','I',size(ifail),'cdiag')
        deallocate(ifail)
      else
#endif
C***************************
C  Serial diagonalisation  *
C***************************
C Serial: call lapack routine chegv/zhegv
C Matrix z used for auxiliary space 

        lwork=5*nh
#ifdef CRAY
        call chegv
#else
        call zhegv
#endif
     .    (1,jobz,uplo,n,h,nh,s,ns,w,z,lwork,fv,info)

C Move eigenvectors to correct storage location
        do i2=1,n
          do i1=1,nh
            z(i1,i2)=h(i1,i2)
          enddo
        enddo
#ifdef MPI
      endif
#endif

      if (info .ne. 0) then
        ierror = 1
        if (Node.eq.0) then
          write(6,*) 'cdiag: ERROR in routine zhegv:'
          if (info.lt.0) then
            write(6,*) -info,'-th argument had an illegal value'
          elseif (info.le.n)   then
#ifdef MPI
            if (TryMemoryIncrease) then
              ierror = -1
              MemoryFactor = 1.5d0*MemoryFactor
              write(6,*) 'Diagonalisation failed to converged'
              write(6,*) 'Memory factor increased to ',MemoryFactor
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
            ierror = -1
            MemoryFactor = 1.5d0*MemoryFactor
          endif
        endif
      endif

c  stop time count
      call timer('cdiag',2)
      end
