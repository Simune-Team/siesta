module m_wallclock
!
! A simple wall-time-stamper
! It should be called only from the master node
!
public :: wallclock

integer, parameter, private :: dp = selected_real_kind(14,200)

integer, private, save     :: last_count
integer, private, save     :: max_count
real(dp), private, save    :: wall_time
integer, private,  save    :: wt
real(dp), private, save    :: count_rate

private

CONTAINS

subroutine wallclock(str)
character(len=*), intent(in) :: str

!
! By using carefully crafted strings, one can get the relevant columns
! to stand out
!
logical :: first = .true.
integer :: count_rate_int
integer        :: count

      if (first) then

         call io_assign(wt)
         open(wt,file='CLOCK',form='formatted',status='unknown')

         CALL system_clock (count_rate=count_rate_int)
         CALL system_clock (count_max=max_count)
         count_rate = real(count_rate_int,kind=dp)
         first = .false.
         CALL system_clock (last_count)
         wall_time = 0.0
         write(unit=wt,fmt="(a,f18.3)") str, wall_time
         RETURN
      endif

CALL system_clock (count)

! Watch out for wrap-around...
! This works if the system counter wrapped around just once...
! ... be liberal in your use of this routine.

if (count < last_count) then
   elapsed_time = (max_count - last_count + count)/count_rate
else
   elapsed_time = (count-last_count)/count_rate
endif

wall_time = wall_time + elapsed_time
write(wt,"(a,f18.3)") str, wall_time
last_count = count

end subroutine wallclock

end module m_wallclock

