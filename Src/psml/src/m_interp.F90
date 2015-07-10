module m_interp

! Default quality parameter for interpolator
integer, public, save                  :: nq = 7

#ifdef __NO_PROC_POINTERS__

! Use a hard-wired interpolator

interface interpolator
   module procedure dpnint1
end interface

public :: interpolator
private :: dpnint1

CONTAINS

#else

! This is the interface that the interpolators
! must have
interface
   subroutine interpolate(nquality,x,y,npts,r,val,debug)

     integer, parameter :: dp = selected_real_kind(10,100)

     integer, intent(in)  :: nquality  ! Quality parameter
     real(dp), intent(in) :: x(*), y(*)
     integer, intent(in)  :: npts    ! Size of x, y arrays
     real(dp), intent(in) :: r
     real(dp), intent(out):: val
     logical, intent(in) :: debug
   end subroutine interpolate
end interface

! 
procedure(interpolate),pointer, public ::  &
                       interpolator => null()
!
! Note that initialization of procedure pointers at declaration
! is a f2008 feature not yet supported by some compilers...
!                       interpolator => dpnint1

public :: set_interpolator, set_default_interpolator
private :: dpnint1

CONTAINS

subroutine set_interpolator(func,nquality)

! Parameter for interpolator's quality
! It might mean different things for different
! interpolators
integer, intent(in) :: nquality

! We should not need to repeat this...
interface
   subroutine func(nquality,x,y,npts,r,val,debug)

     integer, parameter :: dp = selected_real_kind(10,100)

     integer, intent(in)  :: nquality  ! Quality parameter
     real(dp), intent(in) :: x(*), y(*)
     integer, intent(in)  :: npts    ! Size of x, y arrays
     real(dp), intent(in) :: r
     real(dp), intent(out):: val
     logical, intent(in) :: debug
   end subroutine func
end interface

  interpolator => func
  nq = nquality

end subroutine set_interpolator

!
! This routine is needed to work around f2008 issue above
!
subroutine set_default_interpolator()

! Default interpolator and quality parameter
! DRH's dpnint modified by AG, at 7th order
! (Included in this module with permission)
!
  call set_interpolator(dpnint1,7)

end subroutine set_default_interpolator

#endif    /* For systems without procedure pointers */

!
! Copyright (c) 1989-2014 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
! 
! Modified by Alberto Garcia, March 2015
! This routine is included in this module with permission from D.R. Hamann.
!
 subroutine dpnint1(npoly, xx, yy, nn, r, val, debug)

! Modified by Alberto Garcia, March 2015 from routine
! dpnint by D.R. Hamann. 
! Changes:
!   -- A single value is returned
!   -- It can extrapolate, instead of stopping,
!      when called with an abscissa outside the
!      data range.
!   -- If the number of data points is less than
!      npoly+1, npoly is implicitly reduced, without
!      error, and without warning.
!   -- Debug interface 
!
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
 integer ii,imin,imax,iprod,iy,istart,kk,iend

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

!   if (debug) print *, "imin, imax: ", imin, imax

   if(mod(npoly,2)==1) then
    istart=imin-npoly/2
   else if(zz-xx(imin) < xx(imax)-zz) then
     istart=imin-npoly/2
   else
     istart=imax-npoly/2
   end if

   istart = min(istart, nn - npoly)
   istart = max(istart, 1)
   iend = min(istart+npoly,nn)

 !  if (debug) print *, "istart, iend: ", istart, iend
   sum=0.0d0
   do iy=istart,iend
    if(yy(iy)==0.0d0) cycle
    term=yy(iy)
    do iprod=istart, iend
     if(iprod==iy) cycle
     term=term*(zz-xx(iprod))/(xx(iy)-xx(iprod))
    end do
    sum=sum+term
   end do
   val=sum

 end subroutine dpnint1

end module m_interp
