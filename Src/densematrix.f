      module densematrix
C
C  Contains the dense matrix arrays used within SIESTA
C
      use precision

      implicit none

      real(dp), pointer, save :: Haux(:)
      real(dp), pointer, save :: Saux(:)
      real(dp), pointer, save :: psi(:)

      end module densematrix
