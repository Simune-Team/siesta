      module jacobi_davidson
      use matrixop
      use precond,      only : precond_t, setIdPrecond,
     &                         precond_solver_mat,
     &                         precond_solver_vec
      use precision
      use alloc
      use gmres
      implicit none

      PUBLIC :: jd_diag
      PRIVATE

      type    :: jd_opt_t
        real(dp),  pointer :: v0(:) 
        integer            :: kmax
        integer            :: mmin
        integer            :: mmax
        real               :: eps
      end type jd_opt_t

      type(jd_opt_t)       :: jd_opt
      type(precond_t)      :: PC


      contains

      subroutine jd_diag( nuotot, nuo, nspin, neig, fixspin, temp, qtot,
     &                    qs, occtol, numh, listh, H, S, ef, efs,
     &                    Entropy, Dnew, Enew )
      use printMat
      use domain_decom, only : dd_cperm, dd_ncolum

      use parallel,     only : node, nodes

      implicit none
C     Input variables
      integer,   intent(in) :: nuotot, nuo, nspin, neig
      logical,   intent(in) :: fixspin
      real(dp),  intent(in) :: temp, qtot, qs(nspin), occtol
      integer,   intent(in) :: numh(*), listh(*)
      real(dp),  intent(in) :: H(*), S(*)
C     Output variables
      real(dp), intent(out) :: ef, efs(nspin), Entropy
      real(dp), intent(out) :: Dnew(*), Enew(*)
C     Local variables
      integer              :: io, nnz
      integer,     pointer :: ind(:), col(:)
      character*128        :: filename
      type(matrix_t)       :: A, B
      real(dp),    pointer :: w(:), Z(:,:), qo(:)
      real(dp)             :: tt
      logical              :: getD, getPSI

      write(23,*) 'Inside jacobi=', nuotot, nuo, neig
      call pxfflush(23)

      if (Nodes.lt.10) then
        write(filename,'(A7,I1,A3)') 'jacobi.', Node, '.ps'
      else if (Nodes.lt.100) then
        write(filename,'(A7,I2,A3)') 'jacobi.', Node, '.ps'
      else if (Nodes.lt.1000) then
        write(filename,'(A7,I3,A3)') 'jacobi.', Node, '.ps'
      else if (Nodes.lt.10000) then
        write(filename,'(A7,I4,A3)') 'jacobi.', Node, '.ps'
      else
        write(filename,'(A7,I5,A3)') 'jacobi.', Node, '.ps'
      endif

      nullify( ind )
      call re_alloc( ind, 1, nuo+1, 'ind', 'jd_diag' )
      ind(1) = 1
      do io= 1, nuo
        ind(io+1) = ind(io) + numh(io)
      enddo
      nnz = ind(nuo+1)-1

      nullify( col )
      call re_alloc( col, 1, nnz, 'col', 'jd_diag' )
      do io= 1, nnz
        col(io) = dd_cperm(listh(io))
      enddo
      call printMatrix( nuo, nuotot, col, ind, filename, 'MyPart' )

      call create_matrix( nuo, dd_ncolum, nnz, ind, col, H,  A )
      call create_matrix( nuo, dd_ncolum, nnz, ind, col, S,  B )

      nullify( w, Z, qo )
      call re_alloc( w, 1, neig, 'w', 'jd_diag' )
      call re_alloc( Z, 1, nuo, 1, neig, 'Z', 'jd_diag' )
      call re_alloc( qo, 1, neig, 'qo', 'jd_diag' )

C     Set jacobi options
      nullify( jd_opt%v0 )
      call re_alloc( jd_opt%v0, 1, nuo, 'v0', 'jd_diag' )
      tt          = nuotot
      jd_opt%v0   = 1.0/sqrt(tt)
      jd_opt%kmax = neig
      jd_opt%mmax = neig
      jd_opt%mmin = 5
      jd_opt%eps  = 10e-4

C     Set GMRES options and Preconditioner
      gmres_maxiter = nuotot
      gmres_dim     = min(max(257,nuotot/2),nuotot)
      call setIdPrecond( PC, nuo )

      getD   = .true.
      getPSI = .false.

      call jd_diagg( nuo, nspin, neig, nnz, getD, getPSI, fixspin, temp,
     &               qtot, qs, occtol, A, B, ef, efs, Entropy,
     &               w, Z, qo, Dnew, Enew )

      call de_alloc( jd_opt%v0, 'v0', 'jd_diag' )

      call de_alloc( qo, 'qo', 'jd_diag' )
      call de_alloc( Z, 'Z', 'jd_diag' )
      call de_alloc( w, 'w', 'jd_diag' )

      call de_alloc( ind, 'ind', 'jd_diag' )
      call de_alloc( col, 'col', 'jd_diag' )

      end subroutine jd_diag

      subroutine jd_diagg( n, nspin, neig, nnz, getD, getPSI, fixspin,
     &                     temp, qtot, qs, occtol, A, B, ef, efs,
     &                     Entropy, w, Z, qo, Dnew, Enew )
      use m_fermid,      only : fermid, fermispin, stepf
      implicit none
C     Input variables
      integer,        intent(in) :: n, nspin, neig, nnz
      logical,        intent(in) :: getD, getPSI, fixspin
      real(dp)                   :: temp, qtot, qs(nspin), occtol
      type(matrix_t), intent(in) :: A, B
C     Output variables
      real(dp),      intent(out) :: ef, efs(nspin), Entropy
      real(dp),      intent(out), target :: Z(n,neig)
      real(dp),      intent(out) :: w(neig), qo(neig),
     &                              Dnew(nnz), Enew(nnz)
C     Local variables
      integer,         parameter :: nk = 1
      real(dp),        parameter :: wk(1) = (/ 1.0_dp /)
      logical                    :: notfound
      integer                    :: io, ie, j, jo, uneig, nc
      real(dp)                   :: qe, ee, qei, eei
      real(dp),          pointer :: ZZ(:,:)


      call jd_solver( n, neig, A, B, w, Z )
      write(23,*) 'AUTOVALORES de salida', w
      call pxfflush(23)

C Check if we are done
      if (.not.getD) return

C Find new Fermi energy and occupation weights
      if (fixspin) then
        call fermispin( nspin, nspin, nk, wk, neig, neig, 
     .                  w, temp, qs, qo, efs, Entropy )
      else
        call fermid( nspin, nspin, nk, wk, neig, neig, 
     .               w, temp, qtot, qo, ef, Entropy )
      endif

      write(23,*) 'Fermi energy', ef, Entropy
      write(23,*) 'Occupations', qo

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


      write(23,*) 'uneig', uneig
      call pxfflush(23)

      nc = A%nc
#ifdef MPI
      nullify(ZZ)
      call re_alloc( ZZ, 1, nc, 1, uneig, 'ZZ', 'jd_diagg' )
      call reduce_matrix( n, nc, uneig, Z, ZZ )
#else
      ZZ => Z
#endif

      Dnew = 0.0
      Enew = 0.0
      do ie = 1, uneig
        qe = qo(ie)
        ee = qe*w(ie)
        do io = 1, n
          qei = qe*ZZ(io,ie)
          eei = ee*ZZ(io,ie)
          do j = A%ind(io), A%ind(io+1)-1
            jo = A%col(j)
            Dnew(j) = Dnew(j) + qei*ZZ(jo,ie)
            Enew(j) = Enew(j) + eei*ZZ(jo,ie)
          enddo
        enddo
      enddo

      write(23,*) 'Dnew', Dnew
      write(23,*) 'Enew', Enew
      call pxfflush(23)
      call de_alloc( ZZ, 'ZZ', 'jd_diagg' )
      end subroutine jd_diagg

      subroutine jd_solver( n, neig, A, B, lambda, Q )
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      type(matrix_t), intent(in) :: A, B
      integer,        intent(in) :: n, neig
C     Output variables
      real(dp),      intent(out) :: lambda(neig), Q(n,neig)
C     Local variables
      integer              :: kmax, mmin, mmax, it, i, k, m, m2,
     &                        mo, offset, gm_it, gm_it_tot
      real(dp)             :: eps, dotp, norm, tolerance, rteta
      real(dp),    pointer :: t(:), t_tilde(:), V(:,:), Va(:,:),
     &                        Vb(:,:), Vtmp(:,:), Vatmp(:,:),
     &                        Vbtmp(:,:), VAV(:), teta(:), S(:),
     &                        vi(:), via(:), vib(:), vm(:), vma(:),
     &                        vmb(:), s1(:), u(:), p(:),
     &                        ua(:), r(:), si(:), siplus1(:), Z(:,:)

      t   => jd_opt%v0
      kmax = jd_opt%kmax
      mmin = jd_opt%mmin 
      mmax = jd_opt%mmax 
      eps  = jd_opt%eps

      k         = 0
      m         = 0
      gm_it_tot = 0
      nullify( t_tilde, V, Va, Vb, Vtmp, Vatmp, Vbtmp )
      m2 = ((mmax+1)**2+(mmax+1))/2
      call re_alloc( t_tilde, 1, n, 't_tilde', 'jd_solver' )
      call re_alloc( V, 1, n, 1, mmax+1, 'V', 'jd_solver' )
      call re_alloc( Va, 1, n, 1, mmax+1, 'Va', 'jd_solver' )
      call re_alloc( Vb, 1, n, 1, mmax+1, 'Vb', 'jd_solver' )
      call re_alloc( Vtmp, 1, n, 1, mmax+1, 'Vtmp', 'jd_solver' )
      call re_alloc( Vatmp, 1, n, 1, mmax+1, 'Vatmp', 'jd_solver' )
      call re_alloc( Vbtmp, 1, n, 1, mmax+1, 'Vbtmp', 'jd_solver' )
      call re_alloc( VAV, 1, m2, 'VAV', 'jd_solver' )
      call re_alloc( teta, 1, mmax+1, 'teta', 'jd_solver' )
      call re_alloc( S, 1, (mmax+1)**2, 'S', 'jd_solver' )
      call re_alloc( u, 1, n, 'u', 'jd_solver' )
      call re_alloc( p, 1, n, 'p', 'jd_solver' )
      call re_alloc( ua, 1, n, 'ua', 'jd_solver' )
      call re_alloc( r, 1, n, 'r', 'jd_solver' )
      call re_alloc( Z, 1, n, 1, kmax, 'Z', 'jd_solver' )

      ! (1)
      t_tilde = 0.0
      V       = 0.0
      Va      = 0.0
      Vb      = 0.0
      Vtmp    = 0.0
      Vatmp   = 0.0
      Vbtmp   = 0.0
      VAV     = 0.0

      it = 1
      ! (2)
      do while ( k .lt. kmax )
        write(23,*)' '
        write(23,*)'    JACOBI NEW ITER -> m=', m, ' k=', k, ' it=', it
        call pxfflush(23)
        ! (3)
        do i = 1, m
          ! (4)
          vi  =>  V(:,i) 
          vib => Vb(:,i)
          dotp = dvdot( n, vib, t )
          t = t - dotp*vi
        ! (5)  
        enddo

        ! (6) 
        m = m + 1
        call dsmv( n, B, t, t_tilde )
        norm = sqrt(dvdot( n, t_tilde, t ))

        ! (7)
        vm  =>  V(:,m)
        vma => Va(:,m)
        vmb => Vb(:,m)
        vm = t/norm
        call dsmv( n, A, vm, vma )
        vmb = t_tilde/norm

        ! (8) (9) (10)
        offset = (m-1)*m/2
        call ddmtv( n, m, V(1:n,1:m), vma, VAV(offset+1:offset+m) )

        ! (11)
        call compute_eigenpairs( m, VAV, teta, S )
        write(23,*)'      compute_eigenpairs=', teta(1:m)
        call pxfflush(23)

        ! (12)
        s1 => S(1:m)
        call ddmv( n, m, V, s1, u )
        call ddmv( n, m, Vb, s1, p )
        call ddmv( n, m, Va, s1, ua )
        r = ua - teta(1)*p

        ! (13) 
        tolerance = sqrt(dvdot( n, r, r )/dvdot( n, u, u ))
        do while ( tolerance .le. eps )
        write(23,*) '      GOOD EIGENVALUE(', k+1, ')=', teta(1)
        call pxfflush(23)
          ! (14)
          Q(:,k+1)      = u
          Z(:,k+1)      = p
          lambda(k+1)   = teta(1)
          k             = k + 1
          
          ! (15)
          if ( k .eq. kmax ) goto 100

          ! (16)
          mo                  = m
          m                   = m - 1
          VAV(1:(mo*mo+mo)/2) = 0.0

          ! (17)
          if ( m.gt.0 ) then
            do i = 1, m
              ! (18)
              siplus1 => S(1+mo*i:mo*(i+1))
              si      => S(1+m*(i-1):m*i)
              vi      => Vtmp(:,i)
              via     => Vatmp(:,i)
              vib     => Vbtmp(:,i)
              call ddmv( n, mo, V, siplus1, vi )
              call ddmv( n, mo, Va, siplus1, via )
              call ddmv( n, mo, Vb, siplus1, vib )

              ! (19)
              offset = i*(i+1)/2
              VAV(offset) = teta(i+1)
              si          = 0.0
              si(i)       = 1.0
              teta(i)     = teta(i+1)
            ! (20)
            enddo
            V (:,1:m)   = Vtmp (:,1:m)
            Va(:,1:m)   = Vatmp(:,1:m)
            Vb(:,1:m)   = Vbtmp(:,1:m)
          endif
          ! (21)
          u = V(:,1)
          p = Vb(:,1)
          r = Va(:,1) - teta(1)*p

          tolerance = sqrt(dvdot( n, r, r )/dvdot( n, u, u ))

        ! (22)  
        enddo


        ! (23)
        if ( m .gt. mmax ) then
          ! (24)
          VAV = 0.0
          ! (25)
          do i = 2, mmin
            ! (26) 
            si  =>  S(1+m*(i-1):m*i)
            vi  =>  Vtmp(:,i)
            via => Vatmp(:,i)
            vib => Vbtmp(:,i)
            call ddmv( n, m, V, si, vi )
            call ddmv( n, m, Va, si, via )
            call ddmv( n, m, Vb, si, vib )
            offset = i*(i+1)/2
            VAV(offset) = teta(i)
          ! (27)
          enddo
          V (:,2:mmin) =  Vtmp(:,2:mmin)
          Va(:,2:mmin) = Vatmp(:,2:mmin)
          Vb(:,2:mmin) = Vbtmp(:,2:mmin)
          ! (28)
          V(:,1)  = u 
          Va(:,1) = ua
          Vb(:,1) = p
          VAV(1)  = teta(1)
          m       = mmin
        ! (29)
        endif

        ! (30)
        rteta = teta(1) 
        Z(:,k+1) = p
        Q(:,k+1) = u

        ! (31)
        t  = 0.0
        write(23,*)'    PRE jd_correction'
        call pxfflush(23)
        call jd_correction( n, k+1, Z, A, B, rteta, r, gm_it, t )
        write(23,*) '        POS  correction ||t gmres|| = ',
     &              SQRT(dvdot( n, t, t ))
        call pxfflush(23)

        gm_it_tot = gm_it_tot + gm_it
        it = it + 1

      ! (32)  
      end do

100   call de_alloc( Z, 'Z', 'jd_solver' )
      call de_alloc( r, 'r', 'jd_solver' )
      call de_alloc( ua, 'ua', 'jd_solver' )
      call de_alloc( p, 'p', 'jd_solver' )
      call de_alloc( u, 'u', 'jd_solver' )
      call de_alloc( S, 'S', 'jd_solver' )
      call de_alloc( teta, 'teta', 'jd_solver' )
      call de_alloc( VAV, 'VAV', 'jd_solver' )
      call de_alloc( Vbtmp, 'Vbtmp', 'jd_solver' )
      call de_alloc( Vatmp, 'Vatmp', 'jd_solver' )
      call de_alloc( Vtmp, 'Vtmp', 'jd_solver' )
      call de_alloc( Vb, 'Vb', 'jd_solver' )
      call de_alloc( Va, 'Va', 'jd_solver' )
      call de_alloc( V, 'V', 'jd_solver' )
      call de_alloc( t_tilde, 't_tilde', 'jd_solver' )
      end subroutine jd_solver

      subroutine compute_eigenpairs( m, VAV, W, VR )
      use sys,       only : die
      implicit none
C     Input parameters
      integer,    intent(in) :: m
      real(dp)               :: VAV((m*m+m)/2)
C     Output parameters
      real(dp)               :: W(m), VR(m*m)
C     Local Variables
      integer                :: info
      real(dp),      pointer :: A(:), work(:)

      ! LAPACK
      !
      ! DSPEV computes all the eigenvalues and, optionally, eigenvectors
      ! of a complex Hermitian matrix in packed storage
      !
      ! The right eigenvector v(j) of A satisfies
      !                  A * v(j) = lambda(j) * v(j)
      ! where lambda(j) is its eigenvalue.
      ! The left eigenvector u(j) of A satisfies
      !               u(j)**H * A = lambda(j) * u(j)**H
      ! where u(j)**H denotes the conjugate transpose of u(j).
      !
      ! The computed eigenvectors are normalized to have Euclidean norm
      ! equal to 1 and largest component real.
      call re_alloc( A, 1, (m*m+m)/2, 'A', 'compute_eigenpairs' )
      call re_alloc( work, 1, 3*m, 'work', 'compute_eigenpairs' )
      A = VAV

      call DSPEV( 'V', 'U', m, A, W, VR, m, work, info )

      if (info.ne.0) call die( 'ERROR in compute_eigenpairs->DSPEV' )

      call de_alloc( work, 'work', 'compute_eigenpairs' )
      call de_alloc( A, 'A', 'compute_eigenpairs' )

      end subroutine compute_eigenpairs


      subroutine jd_correction( n, ne, Z, A, B, teta, r, gm_it, t )
      use sys,       only : die
      implicit none
C     Input parameters
      integer               :: n, ne
      real(dp)              :: Z(n,ne)
      type(matrix_t)        :: A, B
      real(dp)              :: teta, r(n)
C     Output parameters
      real(dp)              :: t(n)
      integer               :: gm_it
C     Local variables
      integer               :: info
      real(dp),     pointer :: Z_hat(:,:), M(:,:), r_hat(:), gama(:),
     &                         alpha(:), r_tilde(:)
      integer,      pointer :: ipiv(:)

C     Allocate local memory
      call re_alloc( Z_hat, 1, n, 1, ne, 'Z_hat', 'jd_correction' )
      call re_alloc( M, 1, ne, 1, ne, 'M', 'jd_correction' )
      call re_alloc( ipiv, 1, ne, 'ipiv', 'jd_correction' )
      call re_alloc( r_hat, 1, n, 'r_hat', 'jd_correction' )
      call re_alloc( gama, 1, ne, 'gama', 'jd_correction' )
      call re_alloc( alpha, 1, ne, 'alpha', 'jd_correction' )
      call re_alloc( r_tilde, 1, n, 'r_tilde', 'jd_correction' )

      ! (1)
      call precond_solver_mat( n, ne, PC, Z, Z_hat )

      ! (2)
      call ddmtm( n, ne, Z, Z_hat, M )

      ! (3)
C     GETRF computes an LU factorization of a general M-by-N matrix A
C     using partial pivoting with row interchanges.
      call DGETRF( ne, ne, M, ne, ipiv, info )
      if ( info .ne. 0 ) call die( 'ERROR in DGETRF' )


      ! (4)
      ! (5)
      call precond_solver_vec( n, PC, r, r_hat )

      ! (6)
      call ddmtv( n, ne, Z, r_hat, gama )

      ! (7) (8)
C     ZGETRS solves a system of linear equations
C     A * X = B,  A**T * X = B,  or  A**H * X = B
C     with a general N-by-N matrix A using the LU factorization 
C     computed by DGETRF.
      call DGETRS( 'N', ne, 1, M, ne, ipiv, gama, ne, info )
      if (info.ne.0) call die( 'ERROR in jd_correction->DGETRS' )
      alpha = gama

      ! (9) r_tilde = r_hat - Z*alpha
      call ddmv( n, ne, Z, alpha, r_tilde )
C     r_tilde = r_hat - r_tilde
C     we compute -r_tilde
      r_tilde = r_tilde - r_hat
      ! (10)
      write(23,*)'        PRE com_gmres'
      call pxfflush(23)
      call com_gmres( n, ne, A, B, PC, teta, r_tilde, Z, Z_hat,
     &                M, ipiv, t, gm_it )
      write(23,*)'        POS com_gmres'
      call pxfflush(23)

C     Free local memory
      call de_alloc( r_tilde, 'r_tilde', 'jd_correction' )
      call de_alloc( alpha, 'alpha', 'jd_correction' )
      call de_alloc( gama, 'gama', 'jd_correction' )
      call de_alloc( r_hat, 'r_hat', 'jd_correction' )
      call de_alloc( ipiv, 'ipiv', 'jd_correction' )
      call de_alloc( M, 'M', 'jd_correction' )
      call de_alloc( Z_hat, 'Z_hat', 'jd_correction' )

#ifdef KKKK


      ! (10)
      gmres_info%A            = A
      gmres_info%B            = B
      gmres_info%K            = K
      gmres_info%Q_tilde      => Q_tilde    
      gmres_info%Z_tilde_ast  => Z_tilde_ast
      gmres_info%Z_hat        => Z_hat
      gmres_info%M            => M
      gmres_info%teta         = teta
      gmres_info%ne           = ne
      gmres_info%ipiv         => ipiv 

      ierr = 0 
      kryldim = 257
      maxiter = 500
      allocate ( cc(kryldim) )
      allocate ( ss(kryldim) )
      allocate ( rs(kryldim+1) )
      ! TODO: allocate : tamany incorrecte
      allocate ( hh(kryldim + (kryldim * kryldim+1)*2 ) )
      nbnodes = size(t)
      allocate ( kryl(nbnodes,kryldim+1) )
      ! TODO:
      !       - nbnodes ?
      !       - real, complex ?
      !       - treure plot
      allocate ( plot(maxiter) )
      cc   = 0.0
      ss   = 0.0
      rs   = 0.0
      hh   = 0.0
      kryl = 0.0
      plot = 0.0 

#ifdef MPITRACE
        call MPItrace_event ( 1000, 5 )
#endif        
      call gmres_c ( nbnodes,                                           &
                     maxiter, kryldim, eps,                             &
                     -r_tilde,                                          &
                     kryl,                                              &
                     cc, ss, rs, hh,                                    &
                     gmres_info,                                        &
                     t, plot, gmres_iters, ierr )
      if ( ierr.gt.1 ) then
        write(cinfo,*) ierr
        call jd_info  ( 'abnormal termination with ierr = '//cinfo,     &
                        'compute_correction: gmres' )
      endif
#ifdef MPITRACE
        call MPItrace_event ( 1000, 0 )
#endif        

#if DEBUG>=4
     allocate( mxv_wa%y(gmres_info%A%nrows) )
     allocate( mxv_wa%y_hat(gmres_info%A%nrows) )
     allocate( mxv_wa%v_tmp(gmres_info%A%nrows) )
     allocate( mxv_wa%alpha(gmres_info%ne) )
     allocate( mxv_wa%beta (gmres_info%ne) )
     allocate( mxv_wa%gama (gmres_info%ne) )
     call gmres_matrix_dot_vector ( A, B, K, Z_tilde_ast, Z_hat, M, t,  &
                                    teta, v_tmp, mxv_wa, ipiv, Q_tilde )
     call ar_vector_minus_vector ( v_tmp_1, -r_tilde, v_tmp )
     residuo = ar_norm_2(v_tmp_1)/ar_norm_2(-r_tilde)
     write(*,*)'residuo gmres = ', residuo, '; ||v_tmp_1|| = ',         &
                                    ar_norm_2(v_tmp_1)
     deallocate( mxv_wa%y )
     deallocate( mxv_wa%y_hat )
     deallocate( mxv_wa%v_tmp )
     deallocate( mxv_wa%alpha )
     deallocate( mxv_wa%beta  )
     deallocate( mxv_wa%gama )
#endif        


      deallocate ( kryl )               
      deallocate ( cc )
      deallocate ( ss )
      deallocate ( rs )
      deallocate ( hh )
      deallocate ( plot )
#endif
      end subroutine jd_correction













      end module jacobi_davidson
