MODULE siesta_geom
  use precision
  implicit none

  ! Number of atoms in supercell, unit cell
  integer, save :: na_s, na_u

  !unit cell/supercell vectors by columns
  real(dp) :: ucell(3,3), ucell_last(3,3)
  real(dp) :: scell(3,3), scell_last(3,3) 

  ! Atomic coordinates
  real(dp), pointer :: xa(:,:)
  real(dp), pointer :: xalast(:,:)

  ! integer isa(na)           : Species index of each atom
  ! character cisa(na)        : Reference string for each atom
  ! NB cisa is this length in order to contain "siesta:e<isa>"
  ! where isa is the siesta element index, and we allow max 999
  ! such indices 
  integer, pointer  :: isa(:)
  character(len=11), allocatable  :: cisa(:) 

END MODULE siesta_geom
