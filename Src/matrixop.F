      module matrixop
      use precision
      use alloc, only : re_alloc, de_alloc
#ifdef MPI
      use m_mpi_utils,  only : globalize_sum
#endif
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

C
C     DSMV = Double Sparse Matrix Vector
C     y(n) = A(n,n)*x(n)
C
      subroutine dsmv( n, A, x, y )
      implicit none
C     Input variables
      integer,          intent(in) :: n
      type(matrix_t),   intent(in) :: A
      real(dp), target, intent(in) :: x(n)
C     Output variables
      real(dp),       intent(out) :: y(n)
C     Local variables
      integer                     :: ii, jj, kk
      real(dp)                    :: tt
      real(dp),           pointer :: xx(:)
#ifdef MPI
      integer                     :: nc

      nc = A%nc
      call re_alloc( xx, 1, nc, 'xx', 'dsmv' )
      call reduce_array( n, nc, x, xx )
#else
      xx => x
#endif

      do ii= 1, n
        tt = 0.0
        do jj= A%ind(ii), A%ind(ii+1)-1
          kk = A%col(jj)
          tt = tt + A%val(jj)*xx(kk)
        enddo
        y(ii) = tt
      enddo

#ifdef MPI
      call de_alloc( xx, 'xx', 'dsmv' )
#endif
      end subroutine dsmv


C
C     DSMLAMV = Double Sparse Matrix Less Alpha Matrix Vector
C     y(n) = A(n,n)*x(n)
C
      subroutine dsmlamv( n, A, alpha, B, x, y )
      implicit none
C     Input variables
      integer,          intent(in) :: n
      type(matrix_t),   intent(in) :: A, B
      real(dp), target, intent(in) :: alpha, x(n)
C     Output variables
      real(dp),       intent(out) :: y(n)
C     Local variables
      integer                     :: ii, jj, kk
      real(dp)                    :: tt
      real(dp),           pointer :: xx(:)

#ifdef MPI
      integer                     :: nc

      nc = A%nc
      call re_alloc( xx, 1, nc, 'xx', 'dsmv' )
      call reduce_array( n, nc, x, xx )
#else
      xx => x
#endif

      do ii= 1, n
        tt = 0.0
        do jj= A%ind(ii), A%ind(ii+1)-1
          kk = A%col(jj)
          tt = tt + (A%val(jj)-alpha*B%val(jj))*xx(kk)
        enddo
        y(ii) = tt
      enddo

#ifdef MPI
      call de_alloc( xx, 'xx', 'dsmv' )
#endif
      end subroutine dsmlamv

#ifdef MPI
      subroutine reduce_array( n, nc, x, xx )
      use domain_decom, only : dd_ncom, lni, lnb, dd_bsiz, dd_comm
      use parallel,     only : node
      use mpi_siesta
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
      myb => x(lni+1:n)
      ind = n + 1
      do icom= 1, dd_ncom
        PP  = dd_comm(icom)-1
        onb = dd_bsiz(icom)
        if (PP.lt.node) then
          call MPI_Send ( myb, lnb, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, ierror )
          call MPI_Recv ( xx(ind), onb, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, status, ierror )
        else
          call MPI_Recv ( xx(ind), onb, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, status, ierror )
          call MPI_Send ( myb, lnb, MPI_DOUBLE_PRECISION, PP, 0,
     &                    MPI_COMM_WORLD, ierror )
        endif
        ind = ind + onb
      enddo
      end subroutine reduce_array

      subroutine reduce_matrix( n, nc, ne, Z, ZZ )
      use domain_decom, only : dd_ncom, lni, lnb, dd_bsiz, dd_comm
      use parallel,     only : node
      use mpi_siesta
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
      subroutine ddmtv( n, m, A, x, y )
      implicit none
C     Input variables
      integer,    intent(in) :: n, m
      real(dp),   intent(in) :: A(n,m)
      real(dp),   intent(in) :: x(n)
C     Output variables
      real(dp),  intent(out) :: y(m)
C     Local variables
      integer                :: i, j
      real(dp)               :: tt
#ifdef MPI
      real(dp),      pointer :: yt(:)
#endif

      do i= 1, m
        tt = 0.0
        do j= 1, n
          tt = tt + A(j,i)*x(j)
        enddo
        y(i) = tt
      enddo
#ifdef MPI
      call re_alloc( yt, 1, m, 'yt', 'ddmtv' )
      call Globalize_sum( y, yt )
      y = yt
      call de_alloc( yt, 'yt', 'ddmtv' )
#endif
      end subroutine ddmtv

C
C     DDMTM = Double Dense Matrix Transposed Matrix
C     C(nc,nc) = A^t(nr,nc)*B(nr,nc)
C
      subroutine ddmtm( nr, nc, A, B, C )
      implicit none
C     Input variables
      integer,    intent(in) :: nr, nc
      real(dp),   intent(in) :: A(nr,nc), B(nr,nc)
C     Output variables
      real(dp),  intent(out) :: C(nc,nc)
C     Local variables
      integer                :: i, j, k
      real(dp)               :: tt
#ifdef MPI
      real(dp),      pointer :: Ct(:,:)
#endif
      do j= 1, nc
        do i= 1, nc
          tt = 0.0
          do k= 1, nr
            tt = tt + A(k,i)*B(k,j)
          enddo
          C(i,j) = tt
        enddo
      enddo

#ifdef MPI
      call re_alloc( Ct, 1, nc, 1, nc, 'ct', 'ddmtm' )
      call Globalize_sum( C, Ct )
      C = Ct
      call de_alloc( Ct, 'ct', 'ddmtm' )
#endif
      end subroutine ddmtm

C
C     DVDOT = Vector Dot Product
C     r = x(n)*y(n)
C
      function dvdot( n, x, y )
      implicit none
C     Input variables
      integer,  intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
C     Output variables
      real(dp)             :: dvdot
C     Local variables
      integer              :: i
      real(dp)             :: t

      t = 0.0
      do i = 1, n
        t = t + x(i)*y(i)
      enddo

#ifdef MPI
      call Globalize_sum( t, dvdot )
#else
      dvdot = t                    
#endif                          
      return
      end function dvdot

      end module matrixop
