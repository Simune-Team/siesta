module psop_options
!
! Contains options for KB and basis set generation
!
  implicit none

  integer, parameter, private :: dp = selected_real_kind(10,100)

  PUBLIC

  logical :: write_ion_plot_files    ! Write small auxiliary files?
  logical :: new_kb_reference_orbitals ! New scheme for KB reference orbitals
  logical :: debug_kb_generation     ! Write auxiliary files for KB projectors
  logical :: restricted_grid        
  real(dp) :: rmax_ps_check          ! For ps tail checks

end module psop_options
