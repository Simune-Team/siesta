! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module kpoint_scf_m
  !
  ! Contains data structures and routines to deal with the kpoint-grid
  ! for the self-consistent calculation
  ! Other uses (bands, optical, polarization) have their own structures.
  !
  use precision, only : dp
  
  ! The k-point-type
  use kpoint_t_m

  implicit none

  public :: setup_kpoint_scf
  public :: kpoints_scf
  public :: gamma_scf
  
  private

  logical, save :: gamma_scf
  type(kpoint_t), save :: kpoints_scf

contains

  subroutine setup_kpoint_scf( ucell )
    use parallel, only: Node
    use siesta_options, only: writek
    use m_spin, only: TrSym

    real(dp), intent(in) :: ucell(3,3)

    call kpoint_read(kpoints_scf, '', ucell, TrSym)

    gamma_scf = (kpoints_scf%N == 1 .and. &
        dot_product(kpoints_scf%k(:,1),kpoints_scf%k(:,1)) < 1.0e-20_dp)

    ! Quick-return if non-IO
    if ( Node /= 0 ) return

    call kpoint_write_stdout(kpoints_scf, all=writek)
    call kpoint_write_xml(kpoints_scf)
    call kpoint_write_file(kpoints_scf, 'KP')

  end subroutine setup_kpoint_scf
  
end module kpoint_scf_m
