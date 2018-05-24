! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module kpoint_pdos_m
  !
  ! Contains data structures and routines to deal with the kpoint-grid
  ! for the self-consistent calculation
  ! Other uses (bands, optical, polarization) have their own structures.
  !
  use precision, only : dp
  
  ! The k-point-type
  use kpoint_t_m

  implicit none

  public :: setup_kpoint_pdos
  public :: kpoints_pdos
  public :: gamma_pdos

  private

  logical, save :: gamma_pdos
  type(kpoint_t), save :: kpoints_pdos

contains

  subroutine setup_kpoint_pdos( ucell )
    use parallel, only: Node
    use siesta_options, only: writek
    use m_spin, only: TrSym

    real(dp), intent(in) :: ucell(3,3)

    ! First try and read the k-points
    call kpoint_read(kpoints_pdos, 'PDOS', ucell, TrSym)

    if ( kpoints_pdos%method == K_METHOD_NONE ) then
      
      ! The user hasn't specified anything.
      ! This means that we will use the
      call kpoint_delete(kpoints_pdos)
      call kpoint_read(kpoints_pdos, '', ucell, TrSym)

    end if

    gamma_pdos = (kpoints_pdos%N == 1 .and. &
        dot_product(kpoints_pdos%k(:,1),kpoints_pdos%k(:,1)) < 1.0e-20_dp)

    ! Quick-return if non-IO
    if ( Node /= 0 ) return

    ! Write to XML file
    call kpoint_write_stdout(kpoints_pdos, writek, 'PDOS')
    call kpoint_write_xml(kpoints_pdos, 'PDOS')
    call kpoint_write_file(kpoints_pdos, 'PDOS.KP')

  end subroutine setup_kpoint_pdos
  
end module kpoint_pdos_m
