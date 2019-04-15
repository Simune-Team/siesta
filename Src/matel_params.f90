module matel_params

  use precision, only: dp
  
  integer,          parameter :: NQ        =  1024
  integer,          parameter :: NR        =  NQ
  integer,          parameter :: NRTAB     =  1024
  real(dp),         parameter :: Q2CUT     =  2.50e3_dp

  real(dp), parameter :: PI = 3.141592653589_dp
  real(dp), parameter :: QMAX = 2.0_dp * SQRT( Q2CUT )
  real(dp), parameter :: DQ = QMAX / NQ
  real(dp), parameter :: DR = PI / QMAX
  real(dp), parameter :: RMAX = NR * DR
  real(dp), parameter :: DRTAB = RMAX / NRTAB
  
end module matel_params
