! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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
