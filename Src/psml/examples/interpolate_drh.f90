!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! 
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
!
! subroutine dpnint(npoly, xx, yy, nn, r, val, debug)
 subroutine interpolate_drh(npoly, xx, yy, nn, r, val, debug)

! Modified by Alberto Garcia, March 2015 for test purposes

! local polynomial interpolation of data yy on nn points xx
! giving value val on point r
! npoly sets order of polynomial
! xx must be ordered in ascending order
! output interpolated value val on point r

 implicit none

 integer, parameter :: dp=kind(1.0d0)

!Input variables
 real(dp), intent(in) :: xx(*),yy(*)
 real(dp), intent(in) :: r
 real(dp), intent(out) :: val
 integer, intent(in)   ::  nn,npoly
 logical, intent(in)   ::  debug

!Local variables
 real(dp) :: sum,term,zz
 integer ii,imin,imax,iprod,iy,istart,kk

 if(nn<npoly+1) then
   write(6,'(/a,i6,a,i4)') 'dp3int: interpolation error, n=', &
&       nn,'< npoly=',npoly
   stop
 end if

!!$   if(r<xx(1)) then
!!$     write(6,'(/a)') 'dp3int: interpolation error - out of range'
!!$!     stop
!!$   end if
!!$   if(r>xx(nn)) then
!!$     write(6,'(/a)') 'dp3int: interpolation error - out of range'
!!$!     stop
!!$   end if

! interval halving search for xx(ii) points bracketing r

   imin = 1
   imax = nn
   do kk = 1, nn
     ii = (imin + imax) / 2
     if(r>xx(ii)) then
       imin = ii
     else
       imax = ii
     end if
     if(imax - imin .eq. 1) then
       exit
     end if
   end do


   zz=r

   if (debug) print *, "imin, imax: ", imin, imax

   if(mod(npoly,2)==1) then
    istart=imin-npoly/2
   else if(zz-xx(imin) < xx(imax)-zz) then
     istart=imin-npoly/2
   else
     istart=imax-npoly/2
   end if

   istart = min(istart, nn - npoly)
   istart = max(istart, 1)

   if (debug) print *, "istart, iend: ", istart, istart+npoly
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
   val=sum

 end subroutine interpolate_drh

