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

  use m_gauss_fermi, only : G_NF_MIN_kT, G_NF_MAX_kT

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
  integer, parameter :: CC_TYPE_LEN = 17
  integer, parameter :: CC_PART_LEN = 30

! We denote each part by parameters.
! This means that comparing ASCII characters are not an issue
! Furthermore it is clearer in the code for self-explanatory 
! reasons.
  integer, parameter :: CC_PART_EQUI_CIRCLE     =  1
  integer, parameter :: CC_PART_EQUI_LINE       =  2
  integer, parameter :: CC_PART_EQUI_POLES      =  3
  integer, parameter :: CC_PART_L_EQUI_CIRCLE   =  4
  integer, parameter :: CC_PART_L_EQUI_LINE     =  5
  integer, parameter :: CC_PART_L_EQUI_POLES    =  6
  integer, parameter :: CC_PART_R_EQUI_CIRCLE   =  7
  integer, parameter :: CC_PART_R_EQUI_LINE     =  8
  integer, parameter :: CC_PART_R_EQUI_POLES    =  9
  integer, parameter :: CC_PART_NON_EQUI        = 10
  integer, parameter :: CC_PART_TRANSPORT       = 11

! We denote each type of the contour line
! This is merely for book-keeping
  integer, parameter :: CC_TYPE_RES            =   1
  ! The following Fermi Gauss-Quadratures MUST be in success
  ! I.e. NF_<x>kT == 10+x
  integer, parameter :: CC_TYPE_G_NF_MIN       = 4000 ! means G_NF_MIN_kT kT
  integer, parameter :: CC_TYPE_G_NF_MAX       = G_NF_MAX_kT - G_NF_MIN_kT + CC_TYPE_G_NF_MIN ! means G_NF_MAX_kT kT
  integer, parameter :: CC_TYPE_G_NF_0kT       = CC_TYPE_G_NF_MIN - G_NF_MIN_kT ! means 0 kT
  integer, parameter :: CC_TYPE_G_LEGENDRE     = 100
  integer, parameter :: CC_TYPE_G_GEGENBAUER   = 101
  integer, parameter :: CC_TYPE_G_JACOBI       = 102
  integer, parameter :: CC_TYPE_G_CHEBYSHEV_O  = 103
  integer, parameter :: CC_TYPE_G_CHEBYSHEV_C  = 104
  integer, parameter :: CC_TYPE_G_LAGUERRE     = 106
  integer, parameter :: CC_TYPE_G_GEN_LAGUERRE = 107
  integer, parameter :: CC_TYPE_G_HERMITE      = 108
  integer, parameter :: CC_TYPE_SOMMERFELD     = 200
  integer, parameter :: CC_TYPE_SIMP_EXT       = 201
  integer, parameter :: CC_TYPE_SIMP_COMP      = 202
  integer, parameter :: CC_TYPE_SIMP_38        = 203
  integer, parameter :: CC_TYPE_MID            = 204
  integer, parameter :: CC_TYPE_LEFT           = 205
  integer, parameter :: CC_TYPE_RIGHT          = 206


  integer, parameter :: CC_TYPE_TRANSPORT       = 500
! Leave space for the following types...
  integer, parameter :: CC_TYPE_TRANS_PHONON    = 1000

  private :: part2str_int, part2str_ts_ccontour
  interface part2str
     module procedure part2str_int
     module procedure part2str_ts_ccontour
  end interface part2str
  
  private :: type2str_int, type2str_ts_ccontour
  interface type2str
     module procedure type2str_int
     module procedure type2str_ts_ccontour
  end interface type2str

  private :: longtype2str_int, longtype2str_ts_ccontour
  interface longtype2str
     module procedure longtype2str_int
     module procedure longtype2str_ts_ccontour
  end interface longtype2str
  
contains
  
  function part2str_ts_ccontour(c) result(str)
    type(ts_ccontour), intent(in) :: c
    character(len=CC_PART_LEN) :: str
    str = part2str(c%part)
  end function part2str_ts_ccontour
  
  function part2str_int(part) result(str)
    integer, intent(in) :: part
    character(len=CC_PART_LEN) :: str
    if ( part == CC_PART_EQUI_CIRCLE ) then
       str = 'Equilibrium circle'
    else if ( part == CC_PART_EQUI_LINE ) then
       str = 'Equilibrium line'
    else if ( part == CC_PART_EQUI_POLES ) then
       str = 'Equilibrium poles'
    else if ( part == CC_PART_L_EQUI_CIRCLE ) then
       str = 'Left Equilibrium circle'
    else if ( part == CC_PART_L_EQUI_LINE ) then
       str = 'Left Equilibrium line'
    else if ( part == CC_PART_L_EQUI_POLES ) then
       str = 'Left Equilibrium poles'
    else if ( part == CC_PART_R_EQUI_CIRCLE ) then
       str = 'Right Equilibrium circle'
    else if ( part == CC_PART_R_EQUI_LINE ) then
       str = 'Right Equilibrium line'
    else if ( part == CC_PART_R_EQUI_POLES ) then
       str = 'Right Equilibrium poles'
    else if ( part == CC_PART_NON_EQUI ) then
       str = 'Non-equilibrium'
    else if ( part == CC_PART_TRANSPORT ) then
       str = 'Transport'
    else
       call die('Unknown part for the contour')
    end if
  end function part2str_int

  function type2str_ts_ccontour(c) result(str)
    type(ts_ccontour), intent(in) :: c
    character(len=CC_TYPE_LEN) :: str
    str = type2str(c%type)
  end function type2str_ts_ccontour

  function type2str_int(type) result(str)
    integer, intent(in) :: type
    character(len=CC_TYPE_LEN) :: str
    select case ( type )
    case ( CC_TYPE_RES )
       str = 'Residue'
    case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )
       write(str,'(a,i0)') 'G-Fermi_',type-CC_TYPE_G_NF_0kT
    case ( CC_TYPE_G_LEGENDRE )
       str = 'G-Legendre'
    case ( CC_TYPE_G_GEGENBAUER )
       str = 'G-Gegenbauer'
    case ( CC_TYPE_G_JACOBI )
       str = 'G-Jacobi'
    case ( CC_TYPE_G_CHEBYSHEV_O )
       str = 'G-Chebyshev-O'
    case ( CC_TYPE_G_CHEBYSHEV_C )
       str = 'G-Chebyshev-C'
    case ( CC_TYPE_G_LAGUERRE )
       str = 'G-Laguerre'
    case ( CC_TYPE_G_GEN_LAGUERRE )
       str = 'G-Gen.Laguerre'
    case ( CC_TYPE_G_HERMITE )
       str = 'G-Hermite'
    case ( CC_TYPE_SOMMERFELD )
       str = 'Sommerfeld'
    case ( CC_TYPE_SIMP_EXT )
       str = 'Ext. Simpson'
    case ( CC_TYPE_SIMP_COMP )
       str = 'Comp. Simpson'
    case ( CC_TYPE_SIMP_38 )
       str = 'Simpson 3/8'
    case ( CC_TYPE_MID )
       str = 'Mid-rule'
    case ( CC_TYPE_LEFT )
       str = 'Left-rule'
    case ( CC_TYPE_RIGHT )
       str = 'Right-rule'
    case ( CC_TYPE_TRANSPORT )
       str = 'trans'
    case ( CC_TYPE_TRANS_PHONON )
       str = 'phonon'
    case default
       call die('Unknown type for the contour')
    end select
  end function type2str_int

  function longtype2str_ts_ccontour(c) result(str)
    type(ts_ccontour), intent(in) :: c
    character(len=CC_TYPE_LEN*2) :: str
    str = longtype2str(c%type)
  end function longtype2str_ts_ccontour

  function longtype2str_int(type) result(str)
    integer, intent(in) :: type
    character(len=CC_TYPE_LEN*2) :: str
    select case ( type ) 
    case ( CC_TYPE_RES )
       str = 'Residue'
    case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )
       write(str,'(a,'' ('',i0,''kT)'')') 'Gauss-Fermi',type-CC_TYPE_G_NF_0kT
    case ( CC_TYPE_G_LEGENDRE )
       str = 'Gauss-Legendre'
    case ( CC_TYPE_G_GEGENBAUER )
       str = 'Gauss-Gegenbauer'
    case ( CC_TYPE_G_JACOBI )
       str = 'Gauss-Jacobi'
    case ( CC_TYPE_G_CHEBYSHEV_O )
       str = 'Gauss-Chebyshev (open)'
    case ( CC_TYPE_G_CHEBYSHEV_C )
       str = 'Gauss-Chebyshev (closed)'
    case ( CC_TYPE_G_LAGUERRE )
       str = 'Gauss-Laguerre'
    case ( CC_TYPE_G_GEN_LAGUERRE )
       str = 'Gauss-Generalized Laguerre'
    case ( CC_TYPE_G_HERMITE )
       str = 'Gauss-Hermite'
    case ( CC_TYPE_SOMMERFELD )
       str = 'Sommerfeld'
    case ( CC_TYPE_SIMP_EXT )
       str = 'Ext. Simpson'
    case ( CC_TYPE_SIMP_COMP )
       str = 'Comp. Simpson'
    case ( CC_TYPE_SIMP_38 )
       str = 'Simpson 3/8'
    case ( CC_TYPE_MID )
       str = 'Mid-rule'
    case ( CC_TYPE_LEFT )
       str = 'Left-rule'
    case ( CC_TYPE_RIGHT )
       str = 'Right-rule'
    case ( CC_TYPE_TRANSPORT )
       str = 'Transport'
    case ( CC_TYPE_TRANS_PHONON )
       str = 'Phonon-transport'
    case default
       call die('Unknown type for the contour')
    end select
  end function longtype2str_int

end module m_ts_cctype
