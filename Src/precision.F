      module precision

!
!   This module is not really used, but could come handy
!   in the future.
!
!
#ifdef CRAY
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
#else
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
#endif
      end module precision

