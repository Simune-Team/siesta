module m_interpolate
public :: interpolate

CONTAINS

function interpolate(xx,yy,x) result(val)
!
! Interpolate linearly in the (monotonically increasing!) arrays xx and yy
!
integer, parameter :: dp = selected_real_kind(10,100)

real(dp), intent(in) :: xx(:), yy(:)
real(dp), intent(in) :: x
real(dp)             :: val

integer :: i, n

     interface
      subroutine die(str)
      character(len=*), intent(in), optional  :: str
      end subroutine die
     end interface

n = size(xx)
if (size(yy) /= n) call die("Mismatch in array sizes in interpolate")

if ( (x < xx(1)) .or. (x > xx(n))) then
   call die("Interpolate: x not in range")
endif

do i = 2, n
   if (x <= xx(i)) then
      val = yy(i-1) + (x-xx(i-1)) * (yy(i)-yy(i-1))/(xx(i)-xx(i-1))
      exit
   endif
enddo

end function interpolate
end module m_interpolate

#ifdef __TEST__
program test
use m_interpolate, only: interpolate
implicit none
integer, parameter :: dp = selected_real_kind(10,100)
integer, parameter :: n = 10

real(dp) :: xx(n), yy(n)
real(dp) :: x, y
integer  :: i

xx = (/ (0.2*i, i=0,n-1) /)
yy = exp(xx-1)
print "(10f7.3)", xx
print "(10f7.3)", yy

do
   print *, "Enter x:"
   read *, x
   y = interpolate(xx,yy,x)
   print *, x, y
   print *, "Enter y:"
   read *, y
   x = interpolate(yy,xx,y)
   print *, x, y
enddo
end program test
#endif
