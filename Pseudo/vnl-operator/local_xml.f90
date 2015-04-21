module local_xml

  ! Simple module to share the
  ! XML file handle and other data

  use flib_wxml, only: xmlf_t

  integer, parameter :: dp = selected_real_kind(10,100)

  type(xmlf_t), public            :: xf

  real(dp), allocatable :: rl(:), fval(:)
  integer               :: nrl
  real(dp)              :: drl
  logical               :: use_linear_grid

end module local_xml
