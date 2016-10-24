module m_ts_tdir

  implicit none

  ! The transport direction
  ! This corresponds to the cartesian direction
  ! If this is negative, it corresponds to no specific
  ! transport direction.
  integer, save :: ts_tdir = 3
  ! The transport direction
  ! This corresponds to the unit-cell vector index
  integer, save :: ts_tidx = 3

end module m_ts_tdir
