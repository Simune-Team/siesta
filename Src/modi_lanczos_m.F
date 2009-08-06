      module modi_lanczos_m
#ifdef _LANCZOS_
      use precision
      use alloc
      use sys
      use matrixop
      implicit none
      public :: modi_lanczos
      private
      character*1          :: ML_bmat
      character*2          :: ML_which
      integer              :: ML_IDO
      integer              :: ML_N
      integer              :: ML_INFO
      integer              :: ML_iparam(11)
      real(dp)             :: ML_tol
      real(dp),    pointer :: ML_resid(:)
      logical              :: ML_firstime = .true.
      contains
      subroutine modi_lanczos( nuotot, nuo, nspin, fixspin, temp,
     &                         qtot, qs, occtol, numh, listh, H, S,
     &                         neig, ef, efs, Entropy, Dnew, Enew )
      use domain_decom, only : dd_cperm, dd_ncolum
      use m_fermid,     only : fermid, fermispin
      implicit none
C     Input variables
      integer,   intent(in) :: nuotot, nuo, nspin
      logical,   intent(in) :: fixspin
      real(dp),  intent(in) :: temp, qtot, qs(nspin), occtol
      integer,   intent(in) :: numh(*), listh(*)
      real(dp),  intent(in) :: H(*), S(*)
C     Output variables
      integer,  intent(out) :: neig
      real(dp), intent(out) :: ef, efs(nspin), Entropy
      real(dp), intent(out) :: Dnew(*), Enew(*)
C     Local variables
      logical               :: getD, getPSI, notfound
      integer,    parameter :: nk = 1
      real(dp),   parameter :: wk(1) = (/ 1.0_dp /)
      integer               :: ie, io, nc, nnz, ncv, uneig, j, jo
      real(dp)              :: qe, ee, qei, eei
      integer,      pointer :: ind(:), col(:)
      real(dp),     pointer :: eigval(:), eigvec(:,:), qo(:),
     &                         Geigvec(:,:)
      type(matrix_t)        :: A, B

      call timer( 'LANCZOS', 1 )
      getD   = .true.
      getPSI = .false.

      neig = (qtot*1.2)/2.0
      write(23,*) '  Inside modi_lanczos neig=', neig
      write(23,*) '                      qtot=', qtot
      flush(23)

      nullify( ind )
      call re_alloc( ind, 1, nuo+1, 'ind', 'modi_lanczos_m' )
      ind(1) = 1
      do io= 1, nuo
        ind(io+1) = ind(io) + numh(io)
      enddo
      nnz = ind(nuo+1)-1

      nullify( col )
      call re_alloc( col, 1, nnz, 'col', 'modi_lanczos_m' )
      do io= 1, nnz
        col(io) = dd_cperm(listh(io))
      enddo

      call create_matrix( nuo, dd_ncolum, nnz, ind, col, H,  A )
      call create_matrix( nuo, dd_ncolum, nnz, ind, col, S,  B )

      nullify( eigval, eigvec, qo )
      call re_alloc( eigval, 1, neig, 'eigval', 'modi_lanczos_m' )
      call re_alloc( eigvec, 1, nuo, 1, neig, 'eigvec',
     &               'modi_lanczos_m' )
      call re_alloc( qo, 1, neig, 'qo', 'modi_lanczos_m' )

!     Set some parameters for Modified Lanczos
      if (ML_firstime) then
        ML_info      = 0     ! INFO=0 a random vector is generated
        nullify( ML_resid )
        call re_alloc( ML_resid, 1, nuo, 'ML_resid', 'modi_lanczos_m' )
        ML_firstime  = .false.
      else
        ML_info      = 0    ! Use previous residual vector to accelerate method
      endif

      ML_ido       = 0
      ML_bmat      = 'G'     ! Generalized eigenvalue problem
      ML_N         = nuo     ! Dimension of the eigenproblem
      ML_which     = 'SA'    ! SA: smallest algebraic
      ML_tol       = 1.E-4   ! Stopping criterion
      ML_iparam    = 0
      ML_iparam(1) = 1       ! exact shifts with respect to the 
                             ! current Hessenberg matrix
      ML_iparam(3) = 3000    ! maximum number of Arnoldi iterations
      ML_iparam(7) = 2       ! Mode 2 of PSSAUPD: General symmetric

      ML_bmat      = 'I'     ! Generalized eigenvalue problem
      ML_iparam(7) = 1       ! Mode 2 of PSSAUPD: General symmetric

      call mlanczosg( nuotot, nuo, neig, A, B, eigval, eigvec )

C Check if we are done
      if (.not.getD) return

C Find new Fermi energy and occupation weights
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, neig, neig, 
     &                  eigval, temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, nk, wk, neig, neig, 
     &               eigval, temp, qtot, qo, ef, Entropy )
      endif


C     Count the number of usefull eigenvalues/vectors
      uneig    = 1
      notfound = .true.
      do while(notfound)
        if (qo(uneig).gt.occtol) then
          if (uneig+1.gt.neig) then
            notfound = .false.
          else
            uneig = uneig + 1
          endif
        else
          notfound = .false.
          uneig    = uneig - 1
        endif
      enddo

      write(23,*) 'Fermi energy:', ef
      write(23,*) '     Entropy:', Entropy
      write(23,*) ' Occupations:', qo
      write(23,*) '       uneig:', uneig

      nc = A%nc
#ifdef MPI
      nullify(Geigvec)
      call re_alloc( Geigvec, 1, nc, 1, uneig, 'Geigvec',
     &               'modi_lanczos_m' )
      call reduce_matrix( nuo, nc, uneig, eigvec, Geigvec )
#else
      Geigvec => eigvec
#endif

      Dnew(1:nnz) = 0.0
      Enew(1:nnz) = 0.0
      do ie = 1, uneig
        qe = qo(ie)
        ee = qe*eigval(ie)
        do io = 1, nuo
          qei = qe*Geigvec(io,ie)
          eei = ee*Geigvec(io,ie)
          do j = A%ind(io), A%ind(io+1)-1
            jo = A%col(j)
            Dnew(j) = Dnew(j) + qei*Geigvec(jo,ie)
            Enew(j) = Enew(j) + eei*Geigvec(jo,ie)
          enddo
        enddo
      enddo

      write(23,*) '============= eigval ============='
      write(23,*) eigval
      write(23,*) '=================================='

      call de_alloc( Geigvec, 'Geigvec', 'modi_lanczos_m' )

      call re_alloc( qo, 1, neig, 'qo', 'modi_lanczos_m' )
      call de_alloc( eigval, 'eigval', 'modi_lanczos_m' )
      call de_alloc( eigvec, 'eigvec',
     &               'modi_lanczos_m' )

      call de_alloc( col, 'col', 'modi_lanczos_m' )
      call de_alloc( ind, 'ind', 'modi_lanczos_m' )

!      call de_alloc( ML_resid, 'ML_resid', 'modi_lanczos_m' )
      call timer( 'LANCZOS', 2 )
      end subroutine modi_lanczos

      subroutine mlanczosg( n, nuo, nev, A, B, eigval, eigvec )
#ifdef _LANCZOS_
      use SUPERLU
      use lu_factorization, only : LU
#endif
      implicit none
!     Input variables
      integer,        intent(in) :: n, nuo, nev
      type(matrix_t), intent(in) :: A, B
!      type(superlu_t), intent(in) :: LU
!     Output variables
      real(dp),      intent(out) :: eigval(nev), eigvec(nuo,nev)
!     Local  variables
      integer                    :: ncv, lworkl, ipntr(11),
     &                              ierr, nconv, i, j, nr, nc, ldv
      logical                    :: notend, rvec
      real(dp)                   :: sigma, berr
      logical,           pointer :: selct(:)
      real(dp),          pointer :: workd(:), workl(:), V(:,:),
     &                              d(:,:), ax(:), worklu(:), xx(:)
      real(dp)                   :: dvdot
     

      call get_CompRowLoc_Matrix( LU%A, nrow_loc=nr )
      ncv    = min((nev*3)/2,n)
      lworkl = ncv*(ncv+8)
      ldv    = nuo
      nc     = A%nc

!     Allocate memory
      nullify( workd, v, worklu, workl, xx )
      call re_alloc( workd, 1, 3*nuo,   'workd',  'modi_lanczos_m' )
      call re_alloc( V, 1, ldv, 1, ncv, 'V',      'modi_lanczos_m' )
      call re_alloc( workl,  1, lworkl, 'workl',  'modi_lanczos_m' )
      call re_alloc( worklu, 1, nr,     'worklu', 'modi_lanczos_m' )
      call re_alloc( xx,     1, nc,     'xx',     'modi_lanczos_m' )

      notend = .true.
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
      j = 0
      write(23,*) ' '
      write(23,*) '###############################'
      write(23,*) '# STARTING Modified Lanczos Method'
      write(23,*) '# ido           :',  ML_ido
      write(23,*) '# bmat          : ', ML_bmat
      write(23,*) '# N             :',  ML_N
      write(23,*) '# which         : ', ML_which
      write(23,*) '# nev           :',  nev
      write(23,*) '# ncv           :',  ncv
      write(23,*) '# tol           :',  ML_tol
      write(23,*) '# INITIAL iparam:',  ML_iparam
      write(23,*) '###############################'
      write(23,*) ' '
      flush(23)

      berr     = 0.0
      ierr     = 0
      do while(notend)
        j = j + 1
c       %---------------------------------------------%
c       | Repeatedly call the routine PSSAUPD and take| 
c       | actions indicated by parameter IDO until    |
c       | either convergence is indicated or maxitr   |
c       | has been exceeded.                          |
c       %---------------------------------------------%
        call pdsaupd( MPI_COMM_WORLD, ML_IDO, ML_bmat, ML_N, ML_which,
     &                nev, ML_tol, ML_resid, ncv, V, ldv, ML_iparam,
     &                ipntr, workd, workl, lworkl, ML_info )
        if (ML_info.ne.0) then
          write(0,*) ' '
          write(0,*) ' Error in iter = ', j
          write(0,*) ' Error with _saupd, info = ', ML_info,
     &               ' ido = ', ML_ido
          call die( 'Check documentation in _saupd' )
        endif
        if (ML_ido .eq.99) then
          notend = .false.
        else if (ML_ido .eq. -1 .or. ML_ido .eq. 1) then
c         | Perform matrix vector multiplication |
c         |              y <--- A*x              |
c         %--------------------------------------%
          call dsmv( nuo, nc, A, workd(ipntr(1)), xx, workd(ipntr(2)))
        else if (ML_ido .eq. -1 .or. ML_ido .eq. 1) then
c         %--------------------------------------%
c         | Perform  y <--- OP*x = inv[B]*A*x.   |
c         %--------------------------------------%
          call dsmv(nuo, nc, A, workd(ipntr(1)), xx, workd(ipntr(2)))
          call dcopy( nuo, workd(ipntr(2)),  1, workd(ipntr(1)), 1 )

          call reduce_array_lu( 1, nuo, workd(ipntr(2)), nr, worklu )
          call f_pdgssvx( LU%options, LU%A, LU%ScalePermstruct,
     &                    worklu, nr, 1, LU%grid, LU%LUstruct,
     &                    LU%SOLVEstruct, berr, LU%stat, ierr )
          call reduce_array_lu( -1, nr, worklu, nuo, workd(ipntr(2)) )

          if (ierr.ne.0) then
            write(23,*) '  Iter:', j, 'INFO from f_pdgssvx = ', ierr
            flush(23)
            call die( 'Error in LANCZOS' )
          endif

        else if (ML_ido .eq. 2) then
c         %--------------------------------------%
c         | Perform matrix vector multiplication |
c         |              y <--- B*x              |
c         %--------------------------------------%
          call dsmv( nuo, nc, B, workd(ipntr(1)), xx, workd(ipntr(2)))
        else
          write(0,*) ' '
          write(0,*) ' Error with _saupd, info = ', ML_info,
     &               ' ML_ido = ', ML_ido
          write(0,*) ' Check documentation in _saupd '
          write(0,*) ' iparam: ', ML_iparam
          write(0,*) ' '
          call MPI_ABORT( MPI_COMM_WORLD, 0, ierr )
          STOP
        end if 
      enddo

      write(23,*) 'FIRST LOOP DONE ITERS=', j, ML_iparam(3)
      flush(23)

c     %-------------------------------------------%
c     | No fatal errors occurred.                 |
c     | Post-Process using PSSEUPD.               |
c     |                                           |
c     | Computed eigenvalues may be extracted.    |  
c     |                                           |
c     | Eigenvectors may also be computed now if  |
c     | desired.  (indicated by rvec = .true.)    | 
c     %-------------------------------------------%
      rvec = .true.
      nullify( selct, D, ax )
      call re_alloc( selct, 1, ncv, 'selct', 'modi_lanczos_m' )
      call re_alloc( D, 1, ncv, 1, 2, 'D', 'modi_lanczos_m' )
      call re_alloc( ax, 1, nuo, 'ax', 'modi_lanczos_m' )

c     %----------------------------------------------%
c     | Eigenvalues are returned in the first column |
c     | of the two dimensional array D and the       |
c     | corresponding eigenvectors are returned in   |
c     | the first NEV columns of the two dimensional |
c     | array V if requested.  Otherwise, an         |
c     | orthogonal basis for the invariant subspace  |
c     | corresponding to the eigenvalues in D is     |
c     | returned in V.                               |
c     %----------------------------------------------%
      call pdseupd ( MPI_COMM_WORLD, rvec, 'All', selct, D, V, nuo,
     &               sigma, ML_bmat, nuo, ML_which, nev, ML_tol,
     &               ML_resid, ncv, V, nuo, ML_iparam, ipntr, workd,
     &               workl, lworkl, ML_info )
      if ( ML_info.ne.0) then
        write(0,*) ' '
        write(0,*) ' Error with _seupd, info = ', ML_info
        write(0,*) ' Check documentation in _seupd '
        write(0,*) ' iparam: ', ML_iparam
        write(0,*) ' '
        call MPI_ABORT( MPI_COMM_WORLD, 0, ierr )
        STOP
      endif

      nconv = ML_iparam(5)
      do j= 1, nconv
c       %---------------------------%
c       | Compute the residual norm |
c       |                           |
c       |   ||  A*x - lambda*x ||   |
c       |                           |
c       | for the NCONV accurately  |
c       | computed eigenvalues and  |
c       | eigenvectors.  (iparam(5) |
c       | indicates how many are    |
c       | accurate to the requested |
c       | tolerance)                |
c       %---------------------------%
        call dsmv( nuo, nc, A, V(1,j), xx, ax )
        call daxpy( nuo, -D(j,1), V(1,j), 1, ax, 1)
        D(j,2) = sqrt( dvdot( nuo, ax, ax ) )
      enddo
c      %-------------------------------%
c      | Display computed residuals    |
c      %-------------------------------%
C      call pdmout( MPI_COMM_WORLD, 6, nconv, 2, D, ncv, -6,
C     &             'Ritz values and direct residuals' )

      do i= 1, nev
        eigval(i) = D(i,1)
        do j= 1, nuo
          eigvec(j,i) = V(j,i)
        enddo
      enddo
!     Dellocate memory
      call de_alloc( ax,     'ax',     'modi_lanczos_m' )
      call de_alloc( D,      'D',      'modi_lanczos_m' )
      call de_alloc( selct,  'selct',  'modi_lanczos_m' )
      call de_alloc( xx,     'xx',     'modi_lanczos_m' )
      call de_alloc( worklu, 'worklu', 'modi_lanczos_m' )
      call de_alloc( workl,  'workl',  'modi_lanczos_m' )
      call de_alloc( V,      'V',      'modi_lanczos_m' )
      call de_alloc( workd,  'workd',  'modi_lanczos_m' )

!      ML_resid = 0.0
!      do i= 1, nev
!        do j= 1, nuo
!          ML_resid(j) = ML_resid(j) + eigvec(j,i)
!        enddo
!      enddo

!      do j= 1, nuo
!        ML_resid(j) = eigvec(j,1)
!      enddo
      end subroutine mlanczosg

#endif
      end module modi_lanczos_m
