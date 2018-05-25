! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module ts_kpoint_scf_m
  !
  ! Contains data structures and routines to deal with the kpoint-grid
  ! for the self-consistent calculation
  ! Other uses (bands, optical, polarization) have their own structures.
  !
  use precision, only : dp
  
  ! The k-point-type
  use kpoint_t_m

  implicit none

  public :: setup_ts_kpoint_scf
  public :: ts_kpoints_scf
  public :: ts_gamma_scf
  
  private

  logical, save :: ts_gamma_scf
  type(kpoint_t), save :: ts_kpoints_scf

contains

  subroutine setup_ts_kpoint_scf( ucell, kpoints_scf )
    use parallel, only: Node
    use siesta_options, only: writek
    use m_spin, only: TrSym
    use m_ts_global_vars, only : TSmode

    real(dp), intent(in) :: ucell(3,3)
    type(kpoint_t), intent(in) :: kpoints_scf

    call kpoint_read(ts_kpoints_scf, 'TS', ucell, TrSym, process_k_cell=process_k_cell_displ)

    if ( ts_kpoints_scf%method == K_METHOD_NONE ) then

      call kpoint_delete(ts_kpoints_scf)

      if ( TSmode ) then
        ! The user hasn't specified anything.
        ! This means that we will use the default setting from siesta
        call kpoint_read(ts_kpoints_scf, '', ucell, TrSym, process_k_cell=process_k_cell_displ)
      else
        ! To limit memory usage for very high number of k-points
        call kpoint_associate(ts_kpoints_scf, kpoints_scf)
      end if

    end if

    ts_gamma_scf = (ts_kpoints_scf%N == 1 .and. &
        dot_product(ts_kpoints_scf%k(:,1),ts_kpoints_scf%k(:,1)) < 1.0e-20_dp)

    ! Quick-return if non-IO or not a transiesta run
    if ( .not. TSmode ) return
    if ( Node /= 0 ) return

    call kpoint_write_stdout(ts_kpoints_scf, writek, 'transiesta')
    call kpoint_write_xml(ts_kpoints_scf, 'TS')
    call kpoint_write_file(ts_kpoints_scf, 'TS.KP')

  end subroutine setup_ts_kpoint_scf

  subroutine process_k_cell_displ(k_cell, k_displ)
    use m_ts_global_vars, only : TSmode
    use m_ts_tdir, only: ts_tidx

    integer, intent(inout) :: k_cell(3,3)
    real(dp), intent(inout) :: k_displ(3)
    integer :: i

    if ( TSmode .and. ts_tidx > 0 ) then
      i = ts_tidx
      k_cell(:,i) = 0
      k_cell(i,:) = 0
      k_cell(i,i) = 1
      k_displ(i) = 0._dp
    end if

  end subroutine process_k_cell_displ
  
end module ts_kpoint_scf_m
