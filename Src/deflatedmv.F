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

      call medirTiempo( 1, 21 )
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
      call medirTiempo( 1, 22 )
      call dsmlamv( n, nc, A, teta, B, v, vtmp, zz )
      call medirTiempo( 2, 22 )

      ! (12) - (13)   gama = Z*zz
      call medirTiempo( 1, 17 )
      call ddmtv( n, ne, Z, zz, gama, gbuf )
      call medirTiempo( 2, 17 )

      ! (14) (15)   gama = M^(-1)*gama using LU from choleskyFact
      call medirTiempo( 1, 18 )
      call choleskySolve( ne, LU, gama )
      call medirTiempo( 2, 18 )

      ! (16)   zz = zz - Z*gama
      call medirTiempo( 1, 19 )
      call DGEMV( 'N', n, ne, M_ONE, Z, n, gama, 1, ONE, zz, 1 )
      call medirTiempo( 2, 19 )

      call medirTiempo( 2, 21 )
      end subroutine deflatedMV
