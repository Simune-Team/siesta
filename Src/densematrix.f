! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module densematrix
C
C  Contains the dense matrix arrays used within SIESTA
C
      use precision

      implicit none

      real(dp), pointer :: Haux(:)
      real(dp), pointer :: Saux(:)
      real(dp), pointer :: psi(:)

      CONTAINS

      subroutine resetDenseMatrix( )
      use alloc, only : de_alloc
      implicit none

      call de_alloc( Haux, 'Haux', 'densematrix' )
      call de_alloc( Saux, 'Saux', 'densematrix' )
      call de_alloc( psi,  'psi',  'densematrix' )

      end subroutine resetDenseMatrix

      end module densematrix
