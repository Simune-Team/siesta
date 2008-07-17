      MODULE precond
      use precision
      use alloc,     only : re_alloc, de_alloc
      use matrixop,  only : matrix_t
      implicit none
      PUBLIC :: precond_t, setIdPrecond, precond_solver_mat,
     &          precond_solver_vec

      PRIVATE

      type    :: precond_t
        logical                 :: use_diag
        real(dp), pointer       :: ID(:)
        type(matrix_t)          :: L, U
        integer                 :: maxnnzL
        integer                 :: maxnnzU
        integer                 :: fill
        real(dp)                :: tol
        real(dp)                :: alpha
      end type precond_t

      contains

      subroutine setIdPrecond( PREC, n )
      implicit none
C     Input variables
      type(precond_t)     :: PREC
      integer             :: n

      PREC%use_diag = .true.
      nullify(PREC%ID)
      call re_alloc( PREC%ID, 1, n, 'PREC%ID', 'precond' )
      PREC%ID = 1.0
      end subroutine setIdPrecond



      subroutine precond_solver_mat( n, ne, PREC, Z, Z_hat )
      implicit none
C     Input variables
      integer,         intent(in) :: n, ne
      type(precond_t), intent(in) :: PREC
      real(dp),        intent(in) :: Z(n,ne)
C     Output variables
      real(dp),       intent(out) :: Z_hat(n,ne)
C     Local variables
      integer                     :: i, j

      if (PREC%use_diag) then
        do i= 1, ne
          do j= 1, n
            Z_hat(j,i) = Z(j,i)*PREC%ID(j)
          enddo
        enddo
      else

      endif
      end subroutine precond_solver_mat

      subroutine precond_solver_vec( n, PREC, Z, Z_hat )
      implicit none
C     Input variables
      integer,         intent(in) :: n
      type(precond_t), intent(in) :: PREC
      real(dp),        intent(in) :: Z(n)
C     Output variables
      real(dp),       intent(out) :: Z_hat(n)
C     Local variables
      integer                     :: i

      if (PREC%use_diag) then
        do i= 1, n
          Z_hat(i) = Z(i)*PREC%ID(i)
        enddo
      else

      endif
      end subroutine precond_solver_vec

      END MODULE precond
