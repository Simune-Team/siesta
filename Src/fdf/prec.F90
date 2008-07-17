!=====================================================================
! 
! This file is part of the FDF package.
! 
! This module provides precision for integer and reals in FDF library.
! At this moment this module contains precision specification for:
!
!   a) Integer precision     (ip)
!   b) Single Real precision (sp)
!   c) Double Real precision (dp)
!
!
! September 2007
!
!
!=====================================================================

MODULE prec

!
! Precision handling
! Kind parameters 
!
  integer, parameter :: ip = selected_int_kind(9)
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.d0)
END MODULE prec
