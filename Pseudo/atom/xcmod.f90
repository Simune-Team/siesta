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
! This module stores the information about the XC functional to use,
! and provides a function to set it.

! Sample usage:
!   use precision, only: dp
!   use xcmod, only: setxc
!   call setxc( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call atomxc( ... )
!

module xcmod
  use precision, only: dp

  implicit none

  integer, parameter :: maxFunc = 10
  integer,           save :: nXCfunc
  character(len=10), save :: XCauth(MaxFunc), XCfunc(MaxFunc)
  real(dp),          save :: XCweightX(MaxFunc), XCweightC(MaxFunc)

contains

  subroutine setXC( n, func, auth, wx, wc )
    implicit none
    integer,         intent(in):: n       ! Number of functionals
    character(len=*),intent(in):: func(n) ! Functional name labels
    character(len=*),intent(in):: auth(n) ! Functional author labels
    real(dp),        intent(in):: wx(n)   ! Functl weights for exchng
    real(dp),        intent(in):: wc(n)   ! Functl weights for correl
    nXCfunc = n
    XCfunc(1:n) = func(1:n)
    XCauth(1:n) = auth(1:n)
    XCweightX(1:n) = wx(1:n)
    XCweightC(1:n) = wc(1:n)
  end subroutine setXC
end module xcmod
