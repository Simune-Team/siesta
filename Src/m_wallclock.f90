! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_wallclock
!
! A simple wall-time-stamper
! It should be called only from the master node
!

use m_walltime, only: wall_time

public :: wallclock

integer, parameter, private :: dp = selected_real_kind(14,200)

integer, private,  save    :: wt
logical, private,  save    :: first = .true.

private

CONTAINS

subroutine wallclock(str)
character(len=*), intent(in) :: str
!
real(dp) :: t

if (first) then
   call io_assign(wt)
   open(wt,file='CLOCK',form='formatted',status='unknown')
   first = .false.
endif

call wall_time(t)
write(wt,"(a,f18.3)") str, t

end subroutine wallclock

end module m_wallclock

