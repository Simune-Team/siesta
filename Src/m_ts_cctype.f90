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
     sequence
     complex(dp) :: c    ! Contour value
     complex(dp) :: w    ! Contour weight
     integer     :: part ! part of the contour
     integer     :: type ! type of the contour point
  end type ts_ccontour

  ! maximum length of the string that returns the type
  integer, parameter :: CC_TYPE_LEN = 15
  integer, parameter :: CC_PART_LEN = 30

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
  integer, parameter :: CC_TYPE_EQ_RES              =  1
  integer, parameter :: CC_TYPE_EQ_CIRC_G_LEG       =  2
  integer, parameter :: CC_TYPE_EQ_CIRC_G_CH_O      =  3
  integer, parameter :: CC_TYPE_EQ_CIRC_G_CH_C      =  4
  integer, parameter :: CC_TYPE_EQ_FERMI_G_NF       =  5
  integer, parameter :: CC_TYPE_NEQ_SOMMERFELD      =  6
  integer, parameter :: CC_TYPE_NEQ_TAIL_G_NF_0kT   =  7
  integer, parameter :: CC_TYPE_NEQ_TAIL_G_NF_2kT   = 17
  integer, parameter :: CC_TYPE_NEQ_MID_SIMP_EXT    =  8
  integer, parameter :: CC_TYPE_NEQ_MID_SIMP_COMP   =  9
  integer, parameter :: CC_TYPE_NEQ_MID_SIMP_38     = 10
  integer, parameter :: CC_TYPE_NEQ_MID_MID         = 11
  integer, parameter :: CC_TYPE_NEQ_MID_G_LEG       = 12
  integer, parameter :: CC_TYPE_NEQ_MID_G_CH_O      = 13
  integer, parameter :: CC_TYPE_NEQ_MID_G_CH_C      = 14
  integer, parameter :: CC_TYPE_NEQ_TAIL_G_LAGUERRE = 15
  integer, parameter :: CC_TYPE_NEQ_TAIL_G_LEGENDRE = 18
  integer, parameter :: CC_TYPE_NEQ_G_HERMITE       = 16

  ! DEV-notice, always keep TYPE_TRANSPORT as the last element!
  integer, parameter :: CC_TYPE_TRANSPORT       = 18
! Leave space for the following types...
  integer, parameter :: CC_TYPE_TRANS_PHONON    = 100

contains

  function part2str(c) result(str)
    type(ts_ccontour), intent(in) :: c
    character(len=CC_PART_LEN) :: str
    if ( c%part == CC_PART_EQUI ) then
       str = 'Equilibrium'
    else if ( c%part == CC_PART_LEFT_EQUI ) then
       str = 'Left equilibrium'
    else if ( c%part == CC_PART_RIGHT_EQUI ) then
       str = 'Right equilibrium'
    else if ( c%part == CC_PART_NON_EQUI ) then
       str = 'Non-equilibrium'
    else if ( c%part == CC_PART_TRANSPORT ) then
       str = 'Transport'
    end if
  end function part2str

  function type2str(c) result(str)
    type(ts_ccontour), intent(in) :: c
    character(len=CC_TYPE_LEN) :: str
    if ( c%type == CC_TYPE_EQ_RES ) then
       str = 'Eq. Residue'
    else if ( c%type == CC_TYPE_EQ_FERMI_G_NF ) then
       str = 'Eq. G-Fermi'
    else if ( c%type == CC_TYPE_EQ_CIRC_G_LEG ) then
       str = 'Eq.C G-Leg'
    else if ( c%type == CC_TYPE_EQ_CIRC_G_CH_O ) then
       str = 'Eq.C G-Ch-O'
    else if ( c%type == CC_TYPE_EQ_CIRC_G_CH_C ) then
       str = 'Eq.C G-Ch-C'
    else if ( c%type == CC_TYPE_NEQ_SOMMERFELD ) then
       str = 'nEq. SomFeld'
    else if ( c%type == CC_TYPE_NEQ_TAIL_G_NF_0kT ) then
       str = 'nEq. G-Fermi0'
    else if ( c%type == CC_TYPE_NEQ_TAIL_G_NF_2kT ) then
       str = 'nEq. G-Fermi2'
    else if ( c%type == CC_TYPE_NEQ_TAIL_G_LAGUERRE ) then
       str = 'nEq. G-Laguerre'
    else if ( c%type == CC_TYPE_NEQ_G_HERMITE ) then
       str = 'nEq. G-Hermite'
    else if ( c%type == CC_TYPE_NEQ_MID_SIMP_EXT ) then
       str = 'nEq. nF-SimExt'
    else if ( c%type == CC_TYPE_NEQ_MID_SIMP_COMP ) then
       str = 'nEq. nF-SimComp'
    else if ( c%type == CC_TYPE_NEQ_MID_MID ) then
       str = 'nEq. nF-midrule'
    else if ( c%type == CC_TYPE_NEQ_MID_G_LEG ) then
       str = 'nEq. nF-G-Leg'
    else if ( c%type == CC_TYPE_NEQ_MID_G_CH_O ) then
       str = 'nEq. nF-G-Ch-O'
    else if ( c%type == CC_TYPE_NEQ_MID_G_CH_C ) then
       str = 'nEq. nF-G-Ch-C'
    else if ( c%type == CC_TYPE_TRANSPORT ) then
       str = 'trans'
    else if ( c%type == CC_TYPE_TRANS_PHONON ) then
       str = 'phonon'
    end if
  end function type2str

end module m_ts_cctype
