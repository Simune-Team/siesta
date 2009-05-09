      MODULE CG
      use precision
      use alloc,    only : re_alloc, de_alloc
      use precond,  only : precond_t, precond_solver_vec
      use matrixop
      implicit none

      real(dp)            :: cg_eps     = 1.0e-8
      integer             :: cg_dim     = 257
      integer             :: cg_maxiter = 100
      real(dp), parameter :: cg_epsmac  = 1.0d-14


      CONTAINS

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
      logical                     :: conv, FIN2
      integer                     :: i
      real(dp)                    :: vv, alpha, beta, norm, stopcri
      real(dp),           pointer :: rr(:), ww(:), zz(:)

      nullify( rr, ww, zz )
      call re_alloc( rr, 1, n, 'rr', 'com_cg' )
      call re_alloc( ww, 1, n, 'ww', 'com_cg' )
      call re_alloc( zz, 1, n, 'zz', 'com_cg' )
      iters = 0
      rr    = -bb
      ww    = rr
      xx    = 0.0
      vv    = dvdot( n, rr, rr )
      stopcri = cg_eps * sqrt(vv)
      DO i= 1, cg_maxiter
        call cg_MV( n, ne, A, B, PC, teta, ww, Z, Z_hat,
     &              M, ipiv, zz )
        alpha = vv/dvdot( n, ww, zz )
        xx    = xx + alpha*ww
        rr    = rr - alpha*zz
        norm = dvdot( n, rr, rr )

        write(23,*) i, '     CG norm:', sqrt(norm), stopcri
        if( sqrt(norm) .lt. stopcri ) exit

        beta = norm/vv
        ww   = rr + beta*ww
        vv   = norm
      ENDDO
      if (i.gt. cg_maxiter) write(23,*) 'Max number of iterations in CG'


      call de_alloc( zz, 'zz', 'com_cg' )
      call de_alloc( ww, 'ww', 'com_cg' )
      call de_alloc( rr, 'rr', 'com_cg' )
      end subroutine com_cg

      subroutine cg_mv( n, ne, A, B, PC, teta, v, Z, Z_hat,
     &                     M, ipiv, zz )
      use sys,       only : die
      implicit none
C     Input variables
      integer,         intent(in) :: n, ne, ipiv(ne)
      type(matrix_t),  intent(in) :: A, B
      type(precond_t), intent(in) :: PC
      real(dp),        intent(in) :: teta, v(n), Z(n,ne), Z_hat(n,ne),
     &                               M(ne,ne)
C     Output variables
      real(dp),       intent(out) :: zz(n)
C     Local variables
      integer                     :: info
      real(dp),           pointer :: y(:), y_hat(:), gama(:)

      call re_alloc( y, 1, n, 'y', 'cg_mv' )
      call re_alloc( y_hat, 1, n, 'y_hat', 'cg_mv' )
      call re_alloc( gama, 1, ne, 'gama', 'cg_mv' )

      ! (11)
      call dsmlamv( n, A, teta, B, v, y )

      ! (12)
      call precond_solver_vec( n, PC, y, y_hat )

      ! (13)
      call ddmtv( n, ne, Z, y, gama )

      ! (14) (15)
C     DGETRS solves a system of linear equations
C         A * X = B,  A**T * X = B,  or  A**H * X = B
C     with a general N-by-N matrix A using the LU factorization
C     computed by DGETRF.
      call DGETRS ( 'N', ne, 1, M, ne, ipiv, gama, ne, info )
      if (info.ne.0) call die( 'ERROR in cg_mv->DGETRS' )

      ! (16)
      call ddvlmv( n, ne, y_hat, Z_hat, gama, zz )

      call de_alloc( gama, 'gama', 'cg_mv' )
      call de_alloc( y_hat, 'y_hat', 'cg_mv' )
      call de_alloc( y, 'y', 'cg_mv' )

      end subroutine cg_mv

      END MODULE CG
