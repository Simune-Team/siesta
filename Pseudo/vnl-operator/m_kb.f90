module m_kb
  integer, parameter :: dp = selected_real_kind(10,100)
  
type, public :: kb_t
  integer  :: l
  integer  :: seq
  real(dp) :: ekb
  real(dp) :: erefkb
  real(dp) :: dkbcos
  integer  :: nrval
  real(dp), allocatable :: r(:)
  real(dp), allocatable  :: proj(:)
end type kb_t

! Array of structures to hold the KB information
integer, public     :: nprojs
type(kb_t), public, target  :: kbprojs(50)

end module m_kb
