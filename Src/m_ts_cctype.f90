! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_ts_cctype
! Module for containing the complex contour (will also handle the "close"
! to the real axis transport contour).
! The type can easily be streamlined and is containing all relevant information
! in the type.
!
! The important information available is:
!   1. Complex contour point
!   2. Weight of point (for transport contour this is 0.0)
!   3. Which part of the contour
!   4. The type of contour (i.e. which method is used to create it)
! 

  use precision, only : dp

  implicit none

  private :: dp
  public

! Create a type to contain the contour information
  type :: ts_ccontour
     complex(dp) :: c ! Contour value
     complex(dp) :: w ! Contour weight
     integer :: part  ! part of the contour
     integer :: type  ! type of the contour point
  end type ts_ccontour

! We denote each part by parameters.
! This means that comparing ASCII characters are not an issue
! Furthermore it is clearer in the code for self-explanatory 
! reasons.
  integer, parameter :: CC_PART_EQUI        = 1
  integer, parameter :: CC_PART_LEFT_EQUI   = 2
  integer, parameter :: CC_PART_RIGHT_EQUI  = 3
  integer, parameter :: CC_PART_NON_EQUI    = 4
  integer, parameter :: CC_PART_TRANSPORT   = 5

! We denote each type of the contour line
! This is merely for book-keeping
  integer, parameter :: CC_TYPE_RES         = 1
  integer, parameter :: CC_TYPE_FERMI       = 2
  integer, parameter :: CC_TYPE_CIRCLE      = 3
  integer, parameter :: CC_TYPE_NON_EQUI    = 4
  integer, parameter :: CC_TYPE_GAUSS_FERMI = 5
  integer, parameter :: CC_TYPE_GAUSS_QUAD  = 6
  integer, parameter :: CC_TYPE_TRANSPORT   = 7

end module m_ts_cctype
