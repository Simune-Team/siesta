! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
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

      real(dp), pointer, save :: Haux(:)
      real(dp), pointer, save :: Saux(:)
      real(dp), pointer, save :: psi(:)

      CONTAINS

      subroutine resetDenseMatrix( )
      use alloc, only : de_alloc
      implicit none
      CHARACTER(LEN=*), PARAMETER :: MYNAME = 'resetDenseMatrix'

      if (associated(Haux))
     &  call de_alloc( Haux, name='Haux', routine=MYNAME )

      if (associated(Saux))
     &  call de_alloc( Saux, name='Saux', routine=MYNAME )

      if (associated(psi))
     &  call de_alloc( psi, name='psi', routine=MYNAME )

      end subroutine resetDenseMatrix
      end module densematrix
