! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module densematrix
!
!  Contains the dense matrix arrays used within SIESTA
!
  use precision
  
  implicit none

  private
  
  ! Ensure they are initialially nullified
  real(dp), public, pointer :: Haux(:) => null()
  real(dp), public, pointer :: Saux(:) => null()
  real(dp), public, pointer :: psi(:) => null()

  public :: allocDenseMatrix
  public :: resetDenseMatrix

contains

  subroutine allocDenseMatrix(nHaux, nSaux, npsi)
    use alloc, only : re_alloc
    integer, intent(in) :: nHaux, nSaux, npsi

    call re_alloc(Haux, 1, nHaux, 'Haux', 'densematrix', &
        copy=.false., shrink=.false.)
    call re_alloc(Saux, 1, nSaux, 'Saux', 'densematrix', &
        copy=.false., shrink=.false.)
    call re_alloc(psi, 1, npsi, 'psi', 'densematrix', &
        copy=.false., shrink=.false.)

  end subroutine allocDenseMatrix

  subroutine resetDenseMatrix( )
    use alloc, only : de_alloc

    call de_alloc(Haux, 'Haux', 'densematrix')
    call de_alloc(Saux, 'Saux', 'densematrix')
    call de_alloc(psi, 'psi', 'densematrix')
    nullify(Haux, Saux, psi)
    
  end subroutine resetDenseMatrix

end module densematrix
