      module precision
!
!     Precision handling
!
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
!
!     Select precision of certain parts of the program
!     These are set through preprocessor directives. The
!     default behavior is to use single-precision variables
!     (could change in the future)
!
#ifdef BROYDEN_DP
      integer, parameter :: broyden_p = dp
#else
      integer, parameter :: broyden_p = sp
#endif

#ifdef CHARGE_DP
      integer, parameter :: charge_p = dp
#else
      integer, parameter :: charge_p = sp
#endif

      end module precision

