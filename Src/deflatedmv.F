      subroutine deflatedMV( n, ne, kmax, A, B, teta, v, Z,
     &                     LU, zz, work, lwork )
      use precision
      use matrixop
      implicit none
C     Input variables
      integer,         intent(in) :: n, ne, kmax, lwork
      type(matrix_t),  intent(in) :: A, B
      real(dp),        intent(in) :: teta, v(n), Z(n,ne),
     &                               LU((kmax*(kmax+1))/2)
C     Output variables
      real(dp),       intent(out) :: zz(n)
      real(dp),            target :: work(lwork)
C     Local variables
      integer                     :: info, nc, o1, o2
      real(dp),           pointer :: gama(:), gbuf(:),
     &                               vtmp(:)
      real*8,           parameter :: ONE=1.0, M_ONE=-1.0

      nc = A%nc

      o1 = 1
      o2 = o1 + kmax
      gama => work(o1:o2-1)

      o1 = o2
      o2 = o2 + kmax
      gbuf => work(o1:o2-1)

      o1 = o2
      o2 = o2 + nc
      vtmp => work(o1:o2-1)

      ! (11)   zz = (A-teta*B)*v
      call dsmlamv( n, nc, A, teta, B, v, vtmp, zz )

      ! (12) - (13)   gama = Z*zz
      call ddmtv( n, ne, Z, zz, gama, gbuf )

      ! (14) (15)   gama = M^(-1)*gama using LU from choleskyFact
      call choleskySolve( ne, LU, gama )

      ! (16)   zz = zz - Z*gama
      call DGEMV( 'N', n, ne, M_ONE, Z, n, gama, 1, ONE, zz, 1 )
      end subroutine deflatedMV
