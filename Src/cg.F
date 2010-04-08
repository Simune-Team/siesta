      MODULE CG
      use precision
      use alloc,    only : re_alloc, de_alloc
      use precond,  only : precond_t, precond_solver_vec
      use matrixop
      implicit none

      real(dp)            :: cg_epsilon = 1.0e-5
      real(dp)            :: cg_eps     = 1.0e-5
      integer             :: cg_dim     = 257
      integer             :: cg_maxiter = 100
      real(dp), parameter :: cg_epsmac  = 1.0d-14


      CONTAINS

#ifdef KKKK
      subroutine com_cg( n, ne, A, B, PC, teta, bb, Z, Z_hat,
     &                      M, ipiv, xx, iters )
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
      integer                     :: i, nc
      real(dp)                    :: vv, alpha, beta, norm, stopcri
      real(dp),           pointer :: rr(:), ww(:), zz(:), work(:)
      real(dp)                    :: dvdot

      write(23,*) '     Inside CG'
      flush(23)
      nc = A%nc
      nullify( rr, ww, zz, work )
      call re_alloc( rr, 1, n, 'rr', 'com_cg' )
      call re_alloc( ww, 1, n, 'ww', 'com_cg' )
      call re_alloc( zz, 1, n, 'zz', 'com_cg' )
      call re_alloc( work, 1, 2*n+2*ne+nc, 'work', 'minres' )

      iters = 0
      rr    = -bb
      ww    = rr
      xx    = 0.0
      vv    = dvdot( n, rr, rr )
      stopcri = cg_eps * sqrt(vv)
      DO i= 1, cg_maxiter
        write(23,*) i, '     vv:', vv
        flush(23)
        call cg_MV( n, ne, A, B, PC, teta, ww, Z, Z_hat,
     &              M, ipiv, zz, work )
        alpha = vv/dvdot( n, ww, zz )
        xx    = xx + alpha*ww
        rr    = rr - alpha*zz
        norm = dvdot( n, rr, rr )

        write(23,*) i, '     CG norm:', sqrt(norm), stopcri
        flush(23)
        if( sqrt(norm) .lt. stopcri ) exit

        beta = norm/vv
        ww   = rr + beta*ww
        vv   = norm
      ENDDO
      if (i.gt. cg_maxiter) write(23,*) 'Max number of iterations in CG'

      call de_alloc( work, 'work', 'com_cg' )
      call de_alloc( zz,   'zz',   'com_cg' )
      call de_alloc( ww,   'ww',   'com_cg' )
      call de_alloc( rr,   'rr',   'com_cg' )
      end subroutine com_cg
#else
      subroutine com_cg( n, ne, A, B, PC, teta, bb, Z, Z_hat,
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
      integer                     :: i, nc
      real(dp)                    :: a1, a2, alfa, b1, b2, beta,
     &                               stopcri, res
      real(dp),           pointer :: rr(:), ww(:), zz(:), work(:)
      real(dp)                    :: dvdot

      write(23,*) '     Inside CG'
      flush(23)
      nc = A%nc
      nullify( rr, ww, zz, work )
      call re_alloc( rr, 1, n, 'rr', 'com_cg' )
      call re_alloc( ww, 1, n, 'ww', 'com_cg' )
      call re_alloc( zz, 1, n, 'zz', 'com_cg' )
      call re_alloc( work, 1, 2*n+2*ne+nc, 'work', 'minres' )

!     rr = bb - A*xx
      rr = bb

!     ww = -rr
      ww = -rr

!     zz = (A-teta*B)*ww
      call cg_MV( n, ne, A, B, PC, teta, ww, Z, Z_hat,
     &            M, ipiv, zz, work )

!     alfa = (rr*ww)/(ww*zz)
      a1 = dvdot( n, rr, ww )
      a2 = dvdot( n, ww, zz )
      alfa = a1/a2

!     xx = xx + alfa*ww
      xx    = alfa*ww

!     beta = 0.0
      beta = 0.0

!     stopcri = eps*norma(bb)
      res = dvdot( n, bb, bb )
      stopcri = cg_eps * sqrt(res)

      iters = 0
      DO i= 1, cg_maxiter
!       rr = rr - alfa*zz
        rr = rr - alfa*zz

        res = sqrt(dvdot( n, rr, rr ))
        write(23,*) i, '     CG residual norm:', res, stopcri
        flush(23)
        if (res .le. stopcri ) exit

!       beta = (rr*zz)/(ww*zz)
        b1   = dvdot( n, rr, zz )
        b2   = dvdot( n, ww, zz )
        beta = b1/b2

!       ww = -rr +beta*ww
        ww = -rr +beta*ww

!       zz = (A-teta*B)*ww
        call cg_MV( n, ne, A, B, PC, teta, ww, Z, Z_hat,
     &              M, ipiv, zz, work )

!       alfa = (rr*ww)/(ww*zz)
        a1 = dvdot( n, rr, ww )
        a2 = dvdot( n, ww, zz )
        alfa = a1/a2

!       xx = xx + alfa*ww
        xx = xx + alfa*ww
      ENDDO

      if (i.gt. cg_maxiter) write(23,*) 'Max number of iterations in CG'

      call de_alloc( work, 'work', 'com_cg' )
      call de_alloc( zz,   'zz',   'com_cg' )
      call de_alloc( ww,   'ww',   'com_cg' )
      call de_alloc( rr,   'rr',   'com_cg' )
      end subroutine com_cg

#endif
      subroutine cg_MV( n, ne, A, B, PC, teta, v, Z, Z_hat,
     &                  M, ipiv, zz, work )
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
      if (info.ne.0) call die( 'ERROR in cg_MV->DGETRS' )

      ! (16)
      call ddvlmv( n, ne, y_hat, Z_hat, gama, zz )

      end subroutine cg_MV

      END MODULE CG
