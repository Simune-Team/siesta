module m_numbs_neighb
  use precision, only: dp
  implicit none

  integer              :: maxna=200 ! Max. number of neighbor atoms
  integer, allocatable :: jna(:)    ! Index of neighbor atoms
  real(dp), allocatable  :: r2ij(:)    ! Interatomic distances squared
  real(dp), allocatable  :: xij(:,:)   ! Interatomic vectors

end module m_numbs_neighb



