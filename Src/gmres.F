      MODULE GMRES
      use precision
      use alloc,    only : re_alloc, de_alloc
      use precond,  only : precond_t
      use matrixop
      implicit none

      real(dp)            :: gmres_eps     = 1.0e-4
      integer             :: gmres_dim     = 257
      integer             :: gmres_maxiter = 500
      real(dp), parameter :: gmres_epsmac  = 1.0d-14

      real(dp)            :: minres_eps     = 1.0e-4
      integer             :: minres_maxiter = 1000

      CONTAINS

      subroutine com_minres( n, ne, A, B, PC, teta, bb, Z, Z_hat,
     &                       M, ipiv, xx, iters )
      use sys,      only : die
           implicit none
C     Input variables
      integer,         intent(in) :: n, ne, ipiv(ne)
      type(matrix_t),  intent(in) :: A, B
      type(precond_t), intent(in) :: PC
      real(dp),        intent(in) :: teta, bb(n), Z(n,ne), Z_hat(n,ne),
     &                               M(ne,ne)
C     Output variables
      real(dp),       intent(out) :: xx(n)
      integer,        intent(out) :: iters
C     Local variables
      logical                     :: conv
      integer                     :: nc, i, istop
      real(dp)                    :: stopcri, s, oldb, beta, beta1,
     &                               dbar, qnorm, phibar, rhs1, rhs2,
     &                               qrnorm, tnorm2, ynorm2, cs, sn,
     &                               alfa, gamma, delta, epsln, phi,
     &                               gbar, gmin, gmax, oldeps, denom,
     &                               Acond, zz, rtol, epsr, epsx,
     &                               rnorm, epsa, diag, ynorm, anorm
      real(dp),           pointer :: yy(:), r1(:), r2(:), vv(:),
     &                               work(:), ww(:), w1(:), w2(:)
      real(dp),         parameter :: ZERO= 0.0, ONE= 1.0, EPS= 1.0E-14,
     &                               TEN=10.0
      real(dp),          external :: dvdot

      do i= 1, n
        xx(i) = 0.0
      enddo

      nullify( yy, r1 )
      call re_alloc( yy, 1, n, 'yy', 'minres' )
      call re_alloc( r1, 1, n, 'r1', 'minres' )

!     ------------------------------------------------------------------
!     Set up y and v for the first Lanczos vector v1.
!     y  =  beta1 P' v1,  where  P = C**(-1).
!     v is really P' v1.
!     ------------------------------------------------------------------
      call dcopy ( n, bb, 1, yy, 1 )         ! yy = bb
      call dcopy ( n, bb, 1, r1, 1 )         ! r1 = bb

!!!     Ojo si se usa precondicionador....
!!!      if ( precon ) call Msolve( n, b, y )
      beta1  = dvdot( n, bb, yy )

      if (beta1 .lt. ZERO) then    ! M must be indefinite.
        call die( 'ERROR in com_minres-> Precon must be indefinite' )
      end if

      if (beta1 .eq. ZERO) then    ! b = 0 exactly.  Stop with x = 0.
        return
      end if

      beta1  = sqrt( beta1 )       ! Normalize y to get v1 later.

      nc      = A%nc
      nullify( ww, w1, w2, r2, vv, work )
      call re_alloc( ww, 1, n, 'ww', 'minres' )
      call re_alloc( w1, 1, n, 'w1', 'minres' )
      call re_alloc( w2, 1, n, 'w2', 'minres' )
      call re_alloc( r2, 1, n, 'r2', 'minres' )
      call re_alloc( vv, 1, n, 'vv', 'minres' )
      call re_alloc( work, 1, 2*n+2*ne+nc, 'work', 'minres' )
      oldb   = ZERO
      beta   = beta1
      dbar   = ZERO
      epsln  = ZERO
      qrnorm = beta1
      phibar = beta1
      rhs1   = beta1
      rhs2   = ZERO
      tnorm2 = ZERO
      ynorm2 = ZERO
      cs     = - ONE
      sn     = ZERO
      do i= 1, n
        ww(i) = 0.0
      enddo
      do i= 1, n
        w2(i) = 0.0
      enddo
      call dcopy( n, r1, 1, r2, 1 )    ! r2 = r1

      iters = 0
      rtol  = minres_eps
      istop = 0
      DO WHILE (istop.eq.0)
        iters = iters + 1
        s     = 1.0 / beta           ! Normalize previous vector (in y).
        call dcopy( n, yy, 1, vv, 1 )
        call dscal( n, s, vv, 1 )         ! v = vk if P = I

        call gmres_MV( n, ne, A, B, PC, teta, vv, Z, Z_hat,
     &                 M, ipiv, yy, work )
!        call daxpy ( n, (- shift), vv, 1, yy, 1 )
        if (iters .ge. 2) then
          call daxpy ( n, (- beta/oldb), r1, 1, yy, 1 )
        endif

        alfa = dvdot( n, vv, yy )     ! alphak

        call daxpy ( n, (- alfa/beta), r2, 1, yy, 1 )
        call dcopy ( n, r2, 1, r1, 1 )
        call dcopy ( n, yy, 1, r2, 1 )
!        if ( precon ) call Msolve( n, r2, y )

        oldb = beta                     ! oldb = betak
        beta = dvdot( n, r2, yy ) ! beta = betak+1^2
        if (beta .lt. zero) then
           call die( 'ERROR in com_minres-> beta .lt. 0.0' )
        end if

        beta   = sqrt( beta )             ! beta = betak+1
        tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

        if (iters .eq. 1) then              ! Initialize a few things.
          if (beta/beta1 .le. TEN*eps) then ! beta2 = 0 or ~ 0.
            istop = -1                      ! Terminate later.
          end if
          !tnorm2 = alfa**2
           gmax   = abs( alfa )              ! alpha1
           gmin   = gmax                     ! alpha1
        end if

        ! Apply previous rotation Qk-1 to get
        !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
        !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

        oldeps = epsln
        delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
        gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
        epsln  =               sn * beta ! epsln2 = 0         epslnk+1
        dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

        ! Compute the next plane rotation Qk

        gamma  = sqrt( gbar**2 + beta**2 )   ! gammak
        cs     = gbar / gamma                ! ck
        sn     = beta / gamma                ! sk
        phi    = cs * phibar                 ! phik
        phibar = sn * phibar                 ! phibark+1

!       Update  x.

        denom = one/gamma

        do i = 1, n
          w1(i) = w2(i)
          w2(i) = ww(i)
          ww(i) = (vv(i) - oldeps*w1(i) - delta*w2(i))*denom
          xx(i) = xx(i) + phi*ww(i)
        end do

!       Go round again.
        gmax   = max( gmax, gamma )
        gmin   = min( gmin, gamma )
        zz     = rhs1 / gamma
        ynorm2 = zz**2  +  ynorm2
        rhs1   = rhs2  -  delta*zz
        rhs2   =       -  epsln*zz

!       Estimate various norms and test for convergence.
        Anorm  = sqrt( tnorm2 )
        ynorm  = sqrt( ynorm2 )
        epsa   = Anorm * eps
        epsx   = Anorm * ynorm * eps
        epsr   = Anorm * ynorm * rtol
        diag   = gbar
        if (diag .eq. zero) diag = epsa

        qrnorm = phibar
        rnorm  = qrnorm

!       Estimate  cond(A).
!       In this version we look at the diagonals of  R  in the
!       factorization of the lower Hessenberg matrix,  Q * H = R,
!       where H is the tridiagonal matrix from Lanczos with one
!       extra row, beta(k+1) e_k^T.

        Acond  = gmax / gmin

!       See if any of the stopping criteria are satisfied.
!       In rare cases, istop is already -1 from above (Abar = const*I).

        if (istop .eq. 0) then
          if (iters .ge. minres_maxiter ) istop = 5
          if (Acond .ge. 0.1d+0/eps) istop = 4
          if (epsx .ge. beta1     ) istop = 3
          if (qrnorm .le. epsx      ) istop = 2
          if (qrnorm .le. epsr      ) istop = 1
!      write(23,*) 'Loop iters:', iters, ' istop:', istop,
!     &            ' qrnorm:', qrnorm, ' Anorm*ynorm:', Anorm*ynorm
        end if

!       See if it is time to print something.
      ENDDO
      write(23,*) 'FIN com_minres iters:', iters, ' istop:', istop,
     &            ' qrnorm:', qrnorm, ' epsr:', epsr

      call de_alloc( work, 'work', 'minres' )
      call de_alloc( vv, 'vv', 'minres' )
      call de_alloc( r2, 'r2', 'minres' )
      call de_alloc( w2, 'w2', 'minres' )
      call de_alloc( w1, 'w1', 'minres' )
      call de_alloc( ww, 'ww', 'minres' )
      call de_alloc( yy, 'yy', 'minres' )
      call de_alloc( r1, 'r1', 'minres' )
      end subroutine com_minres


      subroutine com_gmres( n, ne, A, B, PC, teta, bb, Z, Z_hat,
     &                      M, ipiv, xx, iters )
      implicit none
C     Input variables
      integer,         intent(in) :: n, ne, ipiv(ne)
      type(matrix_t),  intent(in) :: A, B
      type(precond_t), intent(in) :: PC
      real(dp),        intent(in) :: teta, bb(n), Z(n,ne), Z_hat(n,ne),
     &                               M(ne,ne)
C     Output variables
      real(dp),       intent(out) :: xx(n)
      integer,        intent(out) :: iters
C     Local variables
      logical                     :: conv, FIN2
      integer                     :: ierr, ii, jj, jj1, kk, idx, dim2,
     &                               nc
      real(dp)                    :: raux, gamma, stopcri, ref
      real(dp),           pointer :: kryl(:,:), rs(:), hh(:),
     &                               cc(:), ss(:), work(:)
      real(dp),          external :: dvdot
      real(dp),         parameter :: ONE=1.0, ZERO=0.0

      ierr  = 0
      iters = 0
      conv  = .FALSE.
      xx    = 0.0

!     raux = ||bb||
      raux = SQRT(dvdot( n, bb, bb ))
      ref  = raux

      if (raux .eq. 0.0d0) return

      stopcri = gmres_eps * raux
      dim2    = gmres_dim*(gmres_dim+3)/2
      nc      = A%nc

      nullify( kryl, rs, hh, cc, ss, work )
      call re_alloc( kryl, 1, n, 1, gmres_dim+1, 'kryl', 'com_gmres' )
      call re_alloc( rs, 1, gmres_dim+1, 'rs', 'com_gmres' )
      call re_alloc( hh, 1, dim2, 'hh', 'com_gmres' )
      call re_alloc( cc, 1, gmres_dim, 'cc', 'com_gmres' )
      call re_alloc( ss, 1, gmres_dim, 'ss', 'com_gmres' )
      call re_alloc( work, 1, 2*n+2*gmres_dim+nc, 'work', 'com_gmres' )
!------------------------------------------------- MAIN LOOP
      DO WHILE (.NOT. conv)
        call gmres_MV( n, ne, A, B, PC, teta, xx, Z, Z_hat,
     &                 M, ipiv, kryl, work )

        do ii= 1, n
          kryl(ii,1) = bb(ii) - kryl(ii,1)
        enddo

C       raux = ||kryl(*,1)||
        raux = SQRT(dvdot( n, kryl, kryl))
!        write(23,*) '  GMRES raux:', iters, raux, stopcri

        if (raux.le.stopcri) then
C         The initial guess is the solution
          conv = .TRUE.
        else
!         Initialize 1-st term of the rhs of hessenberg system
          rs(1) = raux
!         Ortonormalize kryl(*,1)
          raux = 1.0/raux
          call DSCAL( N, raux, kryl, 1 )
        endif

        jj   = 0
        idx  = -1
        FIN2 = conv

C       Inner loop. Restarted each kryldim iterations.
        DO WHILE (.NOT. FIN2)
          iters = iters + 1
          jj    = jj + 1
          jj1   = jj + 1
          idx   = idx + jj

          call gmres_MV( n, ne, A, B, PC, teta, kryl(1,jj), Z, Z_hat,
     &                   M, ipiv, kryl(1,jj1), work )
C         Modified Gram-Schmidt
C         For i= 1, j
C           H(i,j) = <v_i, v_j1>
C           v_j1   = v_j1 - H(i,j) * v_i
          call GramSchmidt( n, jj, kryl, kryl(1,jj1), work, hh(idx+1) )

C         H(jj1,jj) = ||kryl(*,jj1)||
          raux = SQRT(dvdot( n, kryl(1,jj1), kryl(1,jj1)))
          hh(idx+jj1) = raux
!          write(23,*) '  GMRES raux:', iters, raux, stopcri

          if (raux.eq.0.0d0) then
            ierr  = 2
            FIN2  = .TRUE.
            conv  = .TRUE.
            idx   = idx - jj
            jj    = jj - 1
            iters = iters - 1
          else
C           Ortonormalize kryl(*,jj1)
            raux = 1.0/raux
            call DSCAL( N, raux, kryl(1,jj1), 1 )

C           Update factorization of H. Perform previous 
C           transformations on jj-th column of H
            do ii= 1, jj-1
              kk   = ii + 1
              raux = hh(idx+ii)

              hh(idx+ii) =  cc(ii)*raux + ss(ii)*hh(idx+kk)
              hh(idx+kk) = -ss(ii)*raux + cc(ii)*hh(idx+kk)
            enddo

            gamma = hh(idx+jj)*hh(idx+jj) + hh(idx+jj1)*hh(idx+jj1)
            gamma = SQRT(gamma)

C           if gamma is zero then take any small
C           value will affect only residual estimate
            if (gamma.eq.0.0) gamma = gmres_epsmac

C           Get next plane rotation
            cc(jj) = hh(idx+jj)/gamma
            ss(jj) = hh(idx+jj1)/gamma

            hh(idx+jj) = cc(jj)*hh(idx+jj) + ss(jj)*hh(idx+jj1)

C           Update the rhs of the LS problem
            rs(jj1) = -ss(jj)*rs(jj)
            rs(jj)  =  cc(jj)*rs(jj)

C           Convergence Test
            raux = ABS( rs(jj1) )
!            write(23,*) '  GMRES raux:', iters, raux, stopcri

            if (raux.le.stopcri) then
!              write(23,*) '         GMRES GOOD...', iters, raux, stopcri
!              flush(23)
              conv = .TRUE.
              FIN2 = .TRUE.
            else
              if (iters.ge.gmres_maxiter) then
!                write(23,*) '         GMRES MAXITER', iters
!                flush(23)
                conv = .TRUE.
                FIN2 = .TRUE.
                ierr = 1
              else if (jj .ge. gmres_dim) then
                FIN2 = .TRUE.
              endif
            endif
          endif
        END DO   ! Krylov Loop End
        if (jj .gt. 0) then
C         Compute y => Solve upper triangular system
          do ii= jj, 2, -1
            rs(ii) = rs(ii)/hh(idx+ii)
            raux   = rs(ii)

            do kk= 1, ii-1
              rs(kk) = rs(kk) - hh(idx+kk) * raux
            enddo

            idx = idx - ii
          enddo
          rs(1) = rs(1) / hh(1)

C         Linear combination of kryl(*,jj)'s to get the solution.
          call DGEMV( 'N', n, jj, ONE, kryl, n, rs, 1, ONE, xx, 1 )
        endif
      END DO 

      call de_alloc( work, 'work', 'com_gmres' )
      call de_alloc( ss,   'ss',   'com_gmres' )
      call de_alloc( cc,   'cc',   'com_gmres' )
      call de_alloc( hh,   'hh',   'com_gmres' )
      call de_alloc( rs,   'rs',   'com_gmres' )
      call de_alloc( kryl, 'kryl', 'com_gmres' )
      end subroutine com_gmres
      END MODULE GMRES

      subroutine gmres_MV( n, ne, A, B, PC, teta, v, Z, Z_hat,
     &                     M, ipiv, zz, work )
      use matrixop
      use sys,      only : die
      use precond,  only : precond_t, precond_solver_vec
      implicit none
C     Input variables
      integer,         intent(in) :: n, ne, ipiv(ne)
      type(matrix_t),  intent(in) :: A, B
      type(precond_t), intent(in) :: PC
      real(dp),        intent(in) :: teta, v(n), Z(n,ne), Z_hat(n,ne),
     &                               M(ne,ne)
C     Output variables
      real(dp),       intent(out) :: zz(n)
      real(dp),            target :: work(*)
C     Local variables
      integer                     :: info, nc
      real(dp),           pointer :: y(:), y_hat(:), gama(:), gbuf(:),
     &                               vv(:)
      call medirTiempo( 1, 9 )

      nc = A%nc

      y     => work(1:n)
      y_hat => work(n+1:2*n)
      gama  => work(2*n+1:2*n+ne)
      gbuf  => work(2*n+ne+1:2*n+2*ne+1)
      vv    => work(2*n+2*ne+1:2*n+2*ne+nc)

      ! (11)
      call dsmlamv( n, nc, A, teta, B, v, vv, y )

      ! (12)
      call precond_solver_vec( n, PC, y, y_hat )

      ! (13)
      call ddmtv( n, ne, Z, y, gama, gbuf )

      ! (14) (15)
C     DGETRS solves a system of linear equations
C         A * X = B,  A**T * X = B,  or  A**H * X = B
C     with a general N-by-N matrix A using the LU factorization
C     computed by DGETRF.
      call DGETRS ( 'N', ne, 1, M, ne, ipiv, gama, ne, info )
      if (info.ne.0) call die( 'ERROR in gmres_mv->DGETRS' )

      ! (16)
      call ddvlmv( n, ne, y_hat, Z_hat, gama, zz )

      call medirTiempo( 2, 9 )

      end subroutine gmres_mv
