!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
subroutine kpoint_convert(ucell,kin,kout,iopt)
! **********************************************************************
! Enables the conversion between Fourier space k-points into reciprocal
! space k-points.
! Created by Nick Papior Andersen, Aug. 2012
! ***************** INPUT **********************************************
! real*8  ucell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! real*8  kin(3)      : k-point in units of [b] or [1/Bohr]
! integer iopt        : -1 => From units of [b] to [1/Bohr]
!                     :  1 => From units of [1/Bohr] to [b]
! ***************** OUTPUT *********************************************
! real*8  kout(3)     : k-point in units of [b] or [1/Bohr]
!
! Allows for conversion between units of reciprocal k-points.
! **********************************************************************
  use precision, only : dp
  use units    , only : Pi
  use sys      , only : die
  
  real(dp), dimension(3,3), intent(in)  :: ucell
  real(dp), dimension(3)  , intent(in)  :: kin
  real(dp), dimension(3)  , intent(out) :: kout
  integer                 , intent(in)  :: iopt
  
! ***********************
! * LOCAL variables     *
! ***********************
  real(dp), dimension(3,3) :: rcell
  integer                  :: i,j
  
  kout(:) = 0.0_dp
  if ( iopt == 1 ) then
     do i=1,3
        do j=1,3
           kout(i) = kout(i) + ucell(j,i)*kin(j)
        end do
        kout(i) = kout(i)/2.0d0/Pi
     end do
  else if ( iopt == -1 ) then
     call reclat(ucell,rcell,0) ! We prefer the consistency of code to not use the ,1 option here
     do i=1,3
        do j=1,3
           kout(i) = kout(i) + rcell(i,j)*kin(j)
        end do
        kout(i) = kout(i)*2.0d0*Pi
     end do
  else
     call die("Wrong option for kpoint_convert! Only 1 or -1 allowed.")
  end if
end subroutine kpoint_convert
