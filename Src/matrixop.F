      module matrixop
      use precision
      use alloc, only : re_alloc, de_alloc
      implicit none

      type    :: matrix_t
        integer           :: nr
        integer           :: nc
        integer           :: nnz
        integer,  pointer :: ind(:)
        integer,  pointer :: col(:)
        real(dp), pointer :: val(:)
      end type matrix_t

      CONTAINS

      subroutine create_matrix( nr, nc, nnz, ind, col, val,  A )
      implicit none
      integer,          intent(in) :: nr, nc, nnz
      integer,  target, intent(in) :: ind(nr+1), col(nnz)
      real(dp), target, intent(in) :: val(nnz)
      type(matrix_t), intent(out) :: A

      A%nr   = nr
      A%nc   = nc
      A%nnz  = nnz
      A%ind => ind
      A%col => col
      A%val => val
      end subroutine create_matrix

      end module matrixop
C
C     DSMV = Double Sparse Matrix Vector
C     y(n) = A(n,n)*x(n)
C
      subroutine dsmv( n, nc, A, x, xx, y )
      use matrixop
      implicit none
C     Input variables
      integer,          intent(in) :: n, nc
      type(matrix_t),   intent(in) :: A
      real(dp),         intent(in) :: x(n)
C     Output variables
      real(dp),       intent(out) :: xx(nc)
      real(dp),       intent(out) :: y(n)
C     Local variables
      integer                     :: ii, jj, kk
      real(dp)                    :: tt

      call medirTiempo( 1, 14 )
#ifdef MPI
      call reduce_array( n, nc, x, xx )
!$OMP PARALLEL DO DEFAULT(NONE)  SCHEDULE(GUIDED)
!$   &  SHARED(n,A,xx,y)
!$   &  PRIVATE(tt,jj,kk)
      do ii= 1, n
        tt = 0.0
        do jj= A%ind(ii), A%ind(ii+1)-1
          kk = A%col(jj)
          tt = tt + A%val(jj)*xx(kk)
        enddo
        y(ii) = tt
      enddo
!$OMP END PARALLEL DO
#else
      do ii= 1, n
        tt = 0.0
        do jj= A%ind(ii), A%ind(ii+1)-1
          kk = A%col(jj)
          tt = tt + A%val(jj)*x(kk)
        enddo
        y(ii) = tt
      enddo

#endif
      call medirTiempo( 2, 14 )
      end subroutine dsmv


C
C     DSMLAMV = Double Sparse Matrix Less Alpha Matrix Vector
C     y(n) = A(n,n)*x(n)
C
      subroutine dsmlamv( n, nc, A, alpha, B, x, xx, y )
      use matrixop
      implicit none
C     Input variables
      integer,          intent(in) :: n, nc
      type(matrix_t),   intent(in) :: A, B
      real(dp), target, intent(in) :: alpha, x(n)
C     Output variables
      real(dp),       intent(out) :: xx(nc)
      real(dp),       intent(out) :: y(n)
C     Local variables
      integer                     :: ii, jj, kk
      real(dp)                    :: tt

      call medirTiempo( 1, 10 )
#ifdef MPI
      call medirTiempo( 1, 11 )
      call reduce_array( n, nc, x, xx )
      call medirTiempo( 2, 11 )
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED)
!$   &  SHARED(n,A,B,xx,y,alpha)
!$   &  PRIVATE(ii,tt,jj,kk)
      do ii= 1, n
        tt = 0.0
        do jj= A%ind(ii), A%ind(ii+1)-1
          kk = A%col(jj)
          tt = tt + (A%val(jj)-alpha*B%val(jj))*xx(kk)
        enddo
        y(ii) = tt
      enddo
#else
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED)
!$   &  SHARED(n,A,B,x,y,alpha)
!$   &  PRIVATE(ii,tt,jj,kk)
      do ii= 1, n
        tt = 0.0
        do jj= A%ind(ii), A%ind(ii+1)-1
          kk = A%col(jj)
          tt = tt + (A%val(jj)-alpha*B%val(jj))*x(kk)
        enddo
        y(ii) = tt
      enddo
#endif
      call medirTiempo( 2, 10 )
      end subroutine dsmlamv

#ifdef MPI
      subroutine reduce_array( n, nc, x, xx )
      use domain_decom, only : dd_ncom, lni, lnb, dd_bsiz, dd_comm
      use parallel,     only : node
      use mpi_siesta
      use precision,    only : dp
      implicit none
C     Input variables
      integer,          intent(in) :: n, nc
      real(dp), target, intent(in) :: x(n)
C     Output variables
      real(dp),        intent(out) :: xx(nc)
C     Local variables
      integer                     :: icom, PP, ind, onb, ierror,
     &                               status(MPI_STATUS_SIZE)
      real(dp),           pointer :: myb(:)

      xx(1:n) = x
      if (dd_ncom.gt.0) then
        myb => x(lni+1:n)
        ind = n + 1
        do icom= 1, dd_ncom
          PP  = dd_comm(icom)-1
          onb = dd_bsiz(icom)
          if (PP.lt.node) then
            call MPI_Send ( myb, lnb, MPI_DOUBLE_PRECISION, PP,
     &                      0, MPI_COMM_WORLD, ierror )
            call MPI_Recv ( xx(ind), onb, MPI_DOUBLE_PRECISION, PP,
     &                      0, MPI_COMM_WORLD, status, ierror )
          else
            call MPI_Recv ( xx(ind), onb, MPI_DOUBLE_PRECISION, PP,
     &                      0, MPI_COMM_WORLD, status, ierror )
            call MPI_Send ( myb, lnb, MPI_DOUBLE_PRECISION, PP,
     &                      0, MPI_COMM_WORLD, ierror )
          endif
          ind = ind + onb
        enddo
      endif
      end subroutine reduce_array

      subroutine reduce_matrix( n, nc, ne, Z, ZZ )
      use domain_decom, only : dd_ncom, lni, lnb, dd_bsiz, dd_comm
      use parallel,     only : node
      use precision,    only : dp
      use mpi_siesta
      use alloc,        only : re_alloc, de_alloc
      implicit none
C     Input variables
      integer,   intent(in) :: n, nc, ne
      real(dp),  intent(in) :: Z(n,ne)
C     Output variables
      real(dp), intent(out) :: ZZ(nc,ne)
C     Local variables
      integer               :: icom, PP, ie, ind, in2, onb, maxbs,
     &                         ierror, status(MPI_STATUS_SIZE)
      real(dp),     pointer :: myb(:,:), buff(:)

      maxbs = 0
      do icom= 1, dd_ncom
        maxbs =max(maxbs,dd_bsiz(icom))
      enddo

      maxbs = maxbs*ne
      nullify(buff,myb)
      call re_alloc( buff, 1, maxbs, 'buff', 'matrixop' )
      call re_alloc( myb, 1, lnb, 1, ne, 'myb', 'matrixop' )

C     Copy the elements that are in the current process
      ZZ(1:n,1:ne) = Z

C     Copy the boundary elements to a contiguos buffer
      myb = Z(lni+1:lni+lnb,1:ne)

      ind = n + 1
      do icom= 1, dd_ncom
        PP  = dd_comm(icom)-1
        onb = dd_bsiz(icom)
        if (PP.lt.node) then
          call MPI_Send ( myb(1,1), lnb*ne, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, ierror )
          call MPI_Recv ( buff, onb*ne, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, status, ierror )
        else
          call MPI_Recv ( buff, onb*ne, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, status, ierror )
          call MPI_Send ( myb(1,1), lnb*ne, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, ierror )
        endif
        
        in2 = 1
        do ie= 1, ne
          ZZ(ind:ind+onb-1,ie) = buff(in2:in2+onb-1)
          in2 = in2 + onb
        enddo
        ind = ind + onb
      enddo
      call de_alloc( myb, 'myb', 'matrixop' )
      call de_alloc( buff, 'buff', 'matrixop' )
      end subroutine reduce_matrix
#endif

C
C     DDMV = Double Dense Matrix Vector
C     y(n) = A(n,m)*x(m)
C
      subroutine ddmv( n, m, A, x, y )
      use matrixop
      use precision,    only : dp
      implicit none
C     Input variables
      integer,    intent(in) :: n, m
      real(dp),   intent(in) :: A(n,m)
      real(dp),   intent(in) :: x(m)
C     Output variables
      real(dp),  intent(out) :: y(n)
C     Local variables
      integer                :: i, j

      y = 0.0
      do i= 1, m
        do j= 1, n
          y(j) = y(j) + A(j,i)*x(i)
        enddo
      enddo
      end subroutine ddmv

C
C     DDVLMV = Double Dense Vector Less Matrix Vector
C     z(n) = x(n) - A(n,m)*y(m)
C
      subroutine ddvlmv( n, m, x, A, y, z )
      use matrixop
      use precision,    only : dp
      implicit none
C     Input variables
      integer,    intent(in) :: n, m
      real(dp),   intent(in) :: x(n)
      real(dp),   intent(in) :: A(n,m)
      real(dp),   intent(in) :: y(m)
C     Output variables
      real(dp),  intent(out) :: z(n)
C     Local variables
      integer                :: i, j

      z = x
      do i= 1, m
        do j= 1, n
          z(j) = z(j) - A(j,i)*y(i)
        enddo
      enddo
      end subroutine ddvlmv


C
C     DDMTV = Double Dense Matrix Transposed Vector
C     y(m) = A^t(n,m)*x(n)
C
      subroutine ddmtv( n, m, A, x, y, yt )
      use matrixop
      use precision,    only : dp
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      integer,    intent(in) :: n, m
      real(dp),   intent(in) :: A(n,m)
      real(dp),   intent(in) :: x(n)
C     Output variables
      real(dp),  intent(out) :: y(m), yt(m)
C     Local variables
      integer                :: i, j
      real(dp)               :: tt
      real(dp),    parameter :: ZERO=0.0, ONE=1.0

#ifdef MPI
      call DGEMV( 'T', n, m, ONE, A, n, x, 1, ZERO, yt, 1 )
      call Globalize_sum( yt, y )
#else
      call DGEMV( 'T', n, m, ONE, A, n, x, 1, ZERO, y, 1 )
#endif
      end subroutine ddmtv

C
C     DDMTM = Double Dense Matrix Transposed Matrix
C     C(nc,nc) = A^t(nr,nc)*B(nr,nc)
C
      subroutine ddmtm( nr, nc, A, B, C, Ct )
      use matrixop
      use precision,    only : dp
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      integer,    intent(in) :: nr, nc
      real(dp),   intent(in) :: A(nr,nc), B(nr,nc)
C     Output variables
      real(dp),  intent(out) :: C(nc,nc), Ct(nc,nc)
C     Local variables
      integer                :: i, j, k
      real(dp)               :: tt

#ifdef MPI
      do j= 1, nc
        do i= 1, nc
          tt = 0.0
          do k= 1, nr
            tt = tt + A(k,i)*B(k,j)
          enddo
          Ct(i,j) = tt
        enddo
      enddo
      call Globalize_sum( Ct, C )
#else
      do j= 1, nc
        do i= 1, nc
          tt = 0.0
          do k= 1, nr
            tt = tt + A(k,i)*B(k,j)
          enddo
          C(i,j) = tt
        enddo
      enddo
#endif
      end subroutine ddmtm

C
C     DVDOT = Vector Dot Product
C     r = x(n)*y(n)
C
      function dvdot( n, x, y )
      use precision,    only : dp
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      integer,  intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
C     Output variables
      real(dp)             :: dvdot
C     Local variables
      integer              :: i
      real(dp)             :: t
      real(dp),   external :: ddot

      t = ddot( n, x, 1, y, 1 )

#ifdef MPI
      call Globalize_sum( t, dvdot )
#else
      dvdot = t                    
#endif                          
      return
      end function dvdot

      subroutine GramSchmidt( N, M, V, t, Vt, gVt )
      use precision,    only : dp
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      integer,     intent(in) :: N, M
      real(dp),    intent(in) :: V(N,M)
C     Output variables
      real(dp)                :: Vt(M), gVt(M)
      real(dp), intent(inout) :: t(N)
C     Local variables
      real(dp),     parameter :: ZERO=0.0, ONE=1.0, MONE=-1.0

      call medirTiempo( 1, 11 )
      call DGEMV( 'T', n, m, ONE, V, n, t, 1, ZERO, Vt, 1 )
#ifdef MPI
      call globalize_sum( Vt, gVt )
      call DGEMV( 'N', n, m, MONE, V, n, gVt, 1, ONE, t, 1 )
#else
      call DGEMV( 'N', n, m, MONE, V, n, Vt, 1, ONE, t, 1 )
#endif
      call medirTiempo( 2, 11 )

      end subroutine GramSchmidt

      subroutine GramSchmidt_B( N, M, VB, V, t, VBt, gVBt )
      use precision,    only : dp
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      integer,     intent(in) :: N, M
      real(dp),    intent(in) :: VB(N,M), V(N,M)
C     Output variables
      real(dp)                :: VBt(M), gVBt(M)
      real(dp), intent(inout) :: t(N)
C     Local variables
      real(dp),     parameter :: ZERO=0.0, ONE=1.0, MONE=-1.0
      real(dp)                :: tin, tout, dvdot, tt
      integer                 :: i

      call medirTiempo( 1, 13 )
#ifdef RGT
      tin = sqrt( dvdot( n, t, t ) )
#endif

#ifdef KKKK
      do i= 1, m
        tt = -dvdot( n, t, Vb(1,i) )
        call DAXPY( n, tt, V(1,i), 1, t, 1 )
      enddo
#else

      call DGEMV( 'T', n, m, ONE, Vb, n, t, 1, ZERO, vbt, 1 )
#ifdef MPI
      call globalize_sum( Vbt, gVbt )
      call DGEMV( 'N', n, m, MONE, V, n, gVbt, 1, ONE, t, 1 )
#else
      call DGEMV( 'N', n, m, MONE, V, n, Vbt, 1, ONE, t, 1 )
#endif
#endif

#ifdef RGT
      tout = sqrt( dvdot( n, t, t ) )
      write(23,*) 'GramSchmidt_B  tout - tin - tout/tin:', tout, tin,
     &            tout/tin
      flush(23)
#endif

      call medirTiempo( 2, 13 )

      end subroutine GramSchmidt_B


      function residual_norm( n, x, y )
      use precision,    only : dp
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
      implicit none
C     Input variables
      integer,  intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
C     Output variables
      real(dp)             :: residual_norm
C     Local variables
      integer              :: i
      real(dp)             :: local(2), global(2)
      real(dp),   external :: ddot

      local(1) = ddot( n, x, 1, x, 1 )
      local(2) = ddot( n, y, 1, y, 1 )
#ifdef MPI
      call globalize_sum( local, global )
      residual_norm = sqrt(global(1)/global(2))
#else
      residual_norm = sqrt(local(1)/local(2))
#endif

#ifdef RGT
      write(23,*) 'residual_norm:', residual_norm
      flush(23)
#endif
      end function residual_norm
