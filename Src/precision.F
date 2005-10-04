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

#ifdef GRID_DP
      integer, parameter :: grid_p = dp
#else
      integer, parameter :: grid_p = sp
#endif

#ifdef ON_DP
      integer, parameter :: on_p = dp
#else
      integer, parameter :: on_p = sp
#endif

      end module precision

