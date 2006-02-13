! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
module units
  ! Define various unit conversion factors from internal units.

  ! internally, siesta works with length: Bohr.
  !                               energy: Rydberg.
  !                                 time: femtosecond

  use precision, only : dp

  implicit none

  real(dp), parameter :: Ang    = 1._dp / 0.529177_dp
  real(dp), parameter :: eV     = 1._dp / 13.60580_dp
  real(dp), parameter :: kBar   = 1._dp / 1.47108e5_dp
  real(dp), parameter :: Kelvin = eV / 11604.45_dp
  real(dp), parameter :: Debye  = 0.393430_dp
  real(dp), parameter :: amu    = 2.133107_dp

! pi to 50 digits
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_dp
  real(dp), parameter :: deg = pi / 180.0_dp

end module units
