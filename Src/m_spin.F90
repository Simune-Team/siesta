! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_spin
  use precision, only: dp
  private

  ! Number of spin components
  integer, save, public           :: nspin

  ! If diagonalization mixes spins, we need to double the block
  ! size for the parallel distribution of matrices. 
  ! In that case, MColl=2,  else MColl=1
  integer, save, public           :: MColl     

  real(dp), pointer, save, public  :: efs(:)
  real(dp), pointer, save, public  :: qs(:)

  public :: init_spin
  public :: fname_spin

CONTAINS

  subroutine init_spin()
    use fdf, only : fdf_get
    use alloc, only: re_alloc
    use parallel, only: IOnode

    implicit none

    logical  :: SPpol, NonCol

    sppol  = fdf_get('SpinPolarized',.false.)
    noncol = fdf_get('NonCollinearSpin',.false.)

!
!   NonCol takes precedence. The printout clarifies
!   whether Up/Down or non-collinear spin is used.
!
    if (NonCol) then

       nspin     = 4
       MColl     = 2
       SPpol     = .false.
    elseif (SPpol) then
       nspin     = 2
       MColl     = 1
    else 
       nspin     = 1
       MColl     = 1
    endif

    if (ionode) then
       write(6,1) 'redata: Non-Collinear-spin run',NonCol
       write(6,1) 'redata: SpinPolarized (Up/Down) run',SPpol
       write(6,2) 'redata: Number of spin components',nspin
    end if

    nullify(efs,qs)
    call re_alloc(efs,1,nspin,name="efs",routine="init_spin")
    call re_alloc(qs,1,nspin,name="qs",routine="init_spin")

1   format(a,t50,'= ',2x,l1)
2   format(a,t50,'= ',2x,i1)


  end subroutine init_spin

  function fname_spin(nspin,ispin,slabel,suffix,basename) result(fname)
    integer, intent(in) :: nspin, ispin
    character(len=*), intent(in), optional :: slabel, suffix, basename
    character(len=200) :: fname
    if ( present(basename) ) then
       if ( nspin == 1 ) then
          fname = trim(basename)
       else
          if ( ispin == 1 ) fname = trim(basename)//'_UP'
          if ( ispin == 2 ) fname = trim(basename)//'_DN'
       end if
    else
       if ( .not. &
            ( present(slabel) .and. present(suffix) ) ) &
            call die('Error in filename input')
       if ( nspin == 1 ) then
          fname = trim(slabel)//'.'//trim(suffix)
       else
          if ( ispin == 1 ) fname = trim(slabel)//'.'//trim(suffix)//'_UP'
          if ( ispin == 2 ) fname = trim(slabel)//'.'//trim(suffix)//'_DN'
       end if
    end if
  end function fname_spin

end module m_spin
