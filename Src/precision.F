! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
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

      public

      end module precision

