! Emulation of the F95 standard cputime, for those compilers
! that do not support it (such as pgf90)
!
! This version for BSD-style system call
!
subroutine cpu_time(t)
real, intent(out)   :: t

real tarray(2)
external etime
real etime

t = etime(tarray)

end subroutine cpu_time

