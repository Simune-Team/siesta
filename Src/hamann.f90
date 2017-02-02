!
! Copyright (c) 1989-2016 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Routines taken from Don Hamann's oncvpsp package

subroutine derivs(mmax, f, a, b, rr, dpr, dppr, dlap)

  ! Computes the first and second derivatives,
  ! and the laplacian 1/r d^2(rf)/dr^2  of a
  ! function f defined on a log grid of the
  ! form r(i) = b * [ exp( a*(i-1) ) - 1]
  !
  ! for which dr/di = a(r+b)
  !
  ! This routine will work also for Hamann-style
  ! grids of the form
  
  ! r(i) = r1 * exp( a*(i-1) ) 
  !
  ! if b is formally set to zero, since dr/di = a*r
  ! in this case (even though r1 (=b) is not zero, obviously)
  
  ! Adapted from a subroutine by Don Hamann
  ! in the oncvpsp package. 
 
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 
 integer, intent(in) :: mmax
 real(dp), intent(in) :: a, b
 real(dp), intent(in) :: rr(mmax)
 real(dp), intent(in) :: f(mmax)
 real(dp), intent(out) :: dpr(mmax), dppr(mmax), dlap(mmax)
 
! local vars
 integer :: i
 real(dp) :: drdi
 real(dp) :: dpn(mmax), dppn(mmax)
 real(dp) :: c11,c12,c13,c14,c15
 real(dp) :: c21,c22,c23,c24,c25
 
 c11 =   2.0d0 / 24.0d0
 c12 = -16.0d0 / 24.0d0
 c13 =   0.0d0 / 24.0d0
 c14 =  16.0d0 / 24.0d0
 c15 =  -2.0d0 / 24.0d0
!
 c21 =   -1.0d0 / 12.0d0
 c22 =   16.0d0 / 12.0d0
 c23 =  -30.0d0 / 12.0d0
 c24 =   16.0d0 / 12.0d0
 c25 =   -1.0d0 / 12.0d0
!
! n derivatives of d
!     
 i=1
 dpn(i) = -25.d0/12.d0*f(i) +4.d0*f(i+1) -3.d0*f(i+2) &
&         +4.d0/3.d0*f(i+3) -1.d0/4.d0*f(i+4)
 dppn(i) = 15.d0/4.d0*f(i) -77.d0/6.d0*f(i+1) +107.d0/6.d0*f(i+2) &
&         -13.d0*f(i+3) +61.d0/12.d0*f(i+4) -5.d0/6.d0*f(i+5)
 i=2
 dpn(i) = -25.d0/12.d0*f(i) +4.d0*f(i+1) -3.d0*f(i+2)  &
&         +4.d0/3.d0*f(i+3) -1.d0/4.d0*f(i+4)
 dppn(i) = 15.d0/4.d0*f(i) -77.d0/6.d0*f(i+1) +107.d0/6.d0*f(i+2) &
&         -13.d0*f(i+3) +61.d0/12.d0*f(i+4) -5.d0/6.d0*f(i+5)
 
 do i = 3, mmax - 2
   dpn(i) =  c11*f(i-2) + c12*f(i-1) + c14*f(i+1) + c15*f(i+2)
   dppn(i) = c21*f(i-2) + c22*f(i-1) + c23*f(i)   + c24*f(i+1) &
&           +c25*f(i+2)
 end do
 
 i=mmax-1
 dpn(i) = +25.d0/12.d0*f(i) -4.d0*f(i-1) +3.d0*f(i-2) &
&         -4.d0/3.d0*f(i-3) +1.d0/4.d0*f(i-4)
 dppn(i) = -15.d0/4.d0*f(i) +77.d0/6.d0*f(i-1) -107.d0/6.d0*f(i-2) &
&          +13.d0*f(i-3) -61.d0/12.d0*f(i-4) +5.d0/6.d0*f(i-5)
 i=mmax
 dpn(i) = +25.d0/12.d0*f(i) -4.d0*f(i-1) +3.d0*f(i-2) &
&         -4.d0/3.d0*f(i-3) +1.d0/4.d0*f(i-4)
 dppn(i) = -15.d0/4.d0*f(i) +77.d0/6.d0*f(i-1) -107.d0/6.d0*f(i-2) &
&          +13.d0*f(i-3) -61.d0/12.d0*f(i-4) +5.d0/6.d0*f(i-5)
 
!
! r derivatives of d
!
 do i = 1, mmax
   drdi = a * (rr(i) + b)
   dpr(i) = dpn(i) / drdi
   dppr(i) = (dppn(i) - a * dpn(i)) / (drdi)**2
   dlap(i) = (dppn(i) + a * dpn(i)) / (drdi)**2
 end do
 
end subroutine derivs


subroutine dpnint(xx, yy, nn, tt, ss, mm)

! local polynomial interpolation of data yy on nn points xx
! giving values ss on mm points tt
! npoly sets order of polynomial
! xx must be ordered in ascending order
! output mm interpolated values ss on points tt

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp) ::  xx(*),yy(*),tt(*),ss(*)
 integer nn,mm

!Local variables
 real(dp) :: sum,term,zz
 integer ii,imin,imax,iprod,iy,istart,jj,kk

! set order of polynomial
 integer, parameter :: npoly=7

 if(nn<npoly+1) then
   write(6,'(/a,i6,a,i4)') 'dpnint: interpolation error, n=', &
&       nn,'< npoly=',npoly
   stop
 end if

! note: output point 1 is skipped in this version because of special
! properties of pp data (ie., log grid)
! this point receives special treatment (see below)

 ss(1:mm)=0.0d0
 imin = 1
 do jj = 2, mm
   if(tt(jj)<xx(1)) then
     write(6,'(/a)') 'dp3int: interpolation error - out of range'
     stop
   end if
   if(tt(jj)>xx(nn)) then
     write(6,'(/a)') 'dpnint: interpolation error - out of range'
     stop
   end if

! interval halving search for xx(ii) points bracketing tt(jj)
   if(jj>2) then
     if(tt(jj)<tt(jj-1)) imin=1
   end if
   imin = 1
   imax = nn
   do kk = 1, nn
     ii = (imin + imax) / 2
     if(tt(jj)>xx(ii)) then
       imin = ii
     else
       imax = ii
     end if
     if(imax - imin .eq. 1) then
       exit
     end if
   end do


   zz=tt(jj)

   if(mod(npoly,2)==1) then
    istart=imin-npoly/2
   else if(zz-xx(imin) < xx(imax)-zz) then
     istart=imin-npoly/2
   else
     istart=imax-npoly/2
   end if

   istart = min(istart, nn - npoly)
   istart = max(istart, 1)

   sum=0.0d0
   do iy=istart,istart+npoly
    if(yy(iy)==0.0d0) cycle
    term=yy(iy)
    do iprod=istart, istart+npoly
     if(iprod==iy) cycle
     term=term*(zz-xx(iprod))/(xx(iy)-xx(iprod))
    end do
    sum=sum+term
   end do
   ss(jj)=sum

! special treatment for the origin
! Do order npoly EXTRAPOLATION to the origin using the points inerpolated
! on the linear grid rather than the log grid, since this represents an
! extrapolation of 1 grid point.  Exrapolation from the innermost points
! of the log grid would represent a huge extrapolation, and round-off
! error would not be acceptable If the fitted function is a polynomial
! of order npoly or less, this is exact.

  if(tt(1)==0.0d0) then

   istart=2
   zz=0.0d0

   sum=0.0d0
   do iy=istart,istart+npoly
    if(ss(iy)==0.0d0) cycle
    term=ss(iy)
    do iprod=istart, istart+npoly
     if(iprod==iy) cycle
     term=term*(zz-tt(iprod))/(tt(iy)-tt(iprod))
    end do
    sum=sum+term
   end do
   ss(1)=sum

  end if

 end do
 return
 end subroutine dpnint
