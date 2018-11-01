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
  public :: resetDenseMatrix_psi

contains

  subroutine allocDenseMatrix(nHaux, nSaux, npsi)
    use alloc, only : re_alloc
    integer, intent(in) :: nHaux, nSaux, npsi

    !> If the arrays are already allocated with the same
    !> bounds nothing will be done
    
    call re_alloc(Haux, 1, nHaux, 'Haux', 'densematrix')
    call re_alloc(Saux, 1, nSaux, 'Saux', 'densematrix')
    call re_alloc(psi, 1, npsi, 'psi', 'densematrix')

  end subroutine allocDenseMatrix

  !> Deallocates auxiliary arrays.
  !> Note that it is safe to call the routine even if
  !> (some) arrays are not associated. Nothing will be
  !> done in that case.
  
  subroutine resetDenseMatrix(keep_psi)
    use alloc, only : de_alloc

    !> This flag is used in connection with the OMM
    !> module: it needs diagon-computed eigenvectors
    !> as seeds for the first few iterations.
    !> [[diagon]] will not deallocate psi in that case
    
    logical, intent(in), optional :: keep_psi

    logical :: keep_psi_allocated

    keep_psi_allocated = .false.
    if (present(keep_psi)) then
       keep_psi_allocated = keep_psi
    endif
       
    call de_alloc(Haux, 'Haux', 'densematrix')
    call de_alloc(Saux, 'Saux', 'densematrix')
    nullify(Haux, Saux)
    
    if (keep_psi_allocated) then
       ! do nothing more
       ! psi will remain allocated
    else
       call de_alloc(psi, 'psi', 'densematrix')
       nullify(psi)
    endif
    
  end subroutine resetDenseMatrix

  subroutine resetDenseMatrix_psi()
    use alloc, only : de_alloc
    
    call de_alloc(psi, 'psi', 'densematrix')
    nullify(psi)
    
  end subroutine resetDenseMatrix_psi

end module densematrix
