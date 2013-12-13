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

  use m_ts_electype, only : NAME_LEN
  use m_ts_io_ctype, only : ts_c_io, CC_METHOD_LEN => c_N
  use precision, only : dp

  use m_gauss_fermi_inf, only : G_NF_MIN_kT, G_NF_MAX_kT

  implicit none

  private :: dp
  public

  type :: ts_ccontour
     sequence
     complex(dp) :: c    ! Contour value
     complex(dp) :: w    ! Contour weight
     integer     :: part ! part of the contour
     integer     :: type ! type of the contour point
  end type ts_ccontour

  ! Create a type to contain the contour information
  type :: ts_eq_c
     type(ts_c_io), pointer :: c_io => null()
     integer :: Ne = 0
     ! the electrode names
     character(len=NAME_LEN), allocatable :: elec(:)
     ! In the same spirit we need to control whether the 
     ! weight has to be altered in any way.
     ! I.e. we have a weight per electrode
     complex(dp), allocatable :: c(:), w(:,:)
  end type ts_eq_c

  ! Create a type to contain the contour information
  type :: ts_neq_c
     type(ts_c_io), pointer :: c_io => null()
     complex(dp), allocatable :: c(:), w(:)
  end type ts_neq_c

  type :: ts_c
     sequence
     logical :: exist = .false.
     complex(dp) :: e
     integer     :: idx(3)
  end type ts_c
  
!  ! maximum length of the string that returns the type
!  integer, parameter :: CC_METHOD_LEN = 17

  ! The following Fermi Gauss-Quadratures MUST be in success
  ! I.e. NF_<x>kT == 10+x
  integer, parameter :: CC_G_NF_MIN         = 4000 ! means G_NF_MIN_kT kT
  integer, parameter :: CC_G_NF_MAX         = G_NF_MAX_kT - G_NF_MIN_kT + CC_G_NF_MIN ! means G_NF_MAX_kT kT
  integer, parameter :: CC_G_NF_0kT         = CC_G_NF_MIN - G_NF_MIN_kT ! means 0 kT
  integer, parameter :: CC_G_LEGENDRE       = 100
  integer, parameter :: CC_TANH_SINH        = 101
  integer, parameter :: CC_SIMP_MIX         = 102
  integer, parameter :: CC_BOOLE_MIX        = 103
  integer, parameter :: CC_MID              = 104


  ! Converts a method to a string format
  interface method2str
     module procedure method2str_int
     module procedure method2str_ts_c_io
  end interface method2str
  private :: method2str_int, method2str_ts_c_io

  ! Converts a method to a string format
  interface longmethod2str
     module procedure longmethod2str_int
     module procedure longmethod2str_ts_c_io
  end interface longmethod2str
  private :: longmethod2str_int, longmethod2str_ts_c_io

  ! Converts a method to the equivalent string format that is required as input
  interface method2input
     module procedure method2input_int
     module procedure method2input_ts_c_io
  end interface method2input
  private :: method2input_int, method2input_ts_c_io

  ! converts a string to the integer method
  interface method
     module procedure method_str
     module procedure method_ts_c_io
  end interface method
  private :: method_str, method_ts_c_io

  interface elec_idx
     module procedure elec_idx_el_eq
     module procedure elec_idx_str_eq
  end interface elec_idx
  private :: elec_idx_el_eq, elec_idx_str_eq

contains

  function elec_idx_el_eq(c,el) result(idx)
    use m_ts_electype
    use fdf, only : leqi
    type(ts_eq_c), intent(in) :: c
    type(Elec), intent(in) :: El
    integer :: idx
    idx = elec_idx(c,name(El))
  end function elec_idx_el_eq

  function elec_idx_str_eq(c,name) result(idx)
    use fdf, only : leqi
    type(ts_eq_c), intent(in) :: c
    character(len=*), intent(in) :: name
    integer :: i, idx
    do i = 1 , size(c%elec)
       if ( leqi(c%elec(i),name) ) then
          idx = i
          return
       end if
    end do
    idx = 0
  end function elec_idx_str_eq

  function method2str_ts_c_io(c) result(str)
    type(ts_c_io), intent(in) :: c
    character(len=CC_METHOD_LEN) :: str
    str = method2str(method(c))
  end function method2str_ts_c_io
  
  function method2str_int(method) result(str)
    integer, intent(in) :: method
    character(len=CC_METHOD_LEN) :: str
    select case ( method )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )
       write(str,'(a,i0)') 'G-Fermi_',method-CC_G_NF_0kT
    case ( CC_G_LEGENDRE )
       str = 'G-Legendre'
    case ( CC_TANH_SINH )
       str = 'Tanh-Sinh'
    case ( CC_SIMP_MIX )
       str = 'Simpson 3/8-3'
    case ( CC_BOOLE_MIX )
       str = 'Boole-Simp-3/8'
    case ( CC_MID )
       str = 'Mid-rule'
    case default
       call die('Unknown method for the contour')
    end select
  end function method2str_int

  function method2input_int(method) result(str)
    integer, intent(in) :: method
    character(len=CC_METHOD_LEN) :: str
    select case ( method )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )
       str = 'g-fermi'
    case ( CC_G_LEGENDRE )
       str = 'g-Legendre'
    case ( CC_TANH_SINH )
       str = 'Tanh-Sinh'
    case ( CC_SIMP_MIX )
       str = 'Simpson-mix'
    case ( CC_BOOLE_MIX )
       str = 'Boole-mix'
    case ( CC_MID )
       str = 'Mid-rule'
    case default
       call die('Unknown method for the contour')
    end select
  end function method2input_int

  function method2input_ts_c_io(c) result(str)
    type(ts_c_io), intent(in) :: c
    character(len=CC_METHOD_LEN) :: str
    str = method2input(method(c))
  end function method2input_ts_c_io

  function longmethod2str_ts_c_io(c) result(str)
    type(ts_c_io), intent(in) :: c
    character(len=CC_METHOD_LEN*2) :: str
    str = longmethod2str(method(c))
  end function longmethod2str_ts_c_io

  function longmethod2str_int(method) result(str)
    integer, intent(in) :: method
    character(len=CC_METHOD_LEN*2) :: str
    select case ( method ) 
    case ( CC_G_NF_MIN:CC_G_NF_MAX )
       write(str,'(a,'' ('',i0,''kT)'')') 'Gauss-Fermi',method-CC_G_NF_0kT
    case ( CC_G_LEGENDRE )
       str = 'Gauss-Legendre'
    case ( CC_TANH_SINH )
       str = 'Tanh-Sinh'
    case ( CC_SIMP_MIX )
       str = 'Simpson-mix'
    case ( CC_BOOLE_MIX )
       str = 'Boole-mix'
    case ( CC_MID )
       str = 'Mid-rule'
    case default
       call die('Unknown method for the contour')
    end select
  end function longmethod2str_int

  function method_str(str) result(method)
    use fdf, only : leqi
    character(len=*), intent(in) :: str
    integer :: method
    character(len=20) :: tmp
    integer :: i
    if ( leqi(str,'g-legendre') ) then
       method = CC_G_LEGENDRE
    else if ( leqi(str,'tanh-sinh') ) then
       method = CC_TANH_SINH
    else if ( leqi(str,'simpson-mix') ) then
       method = CC_SIMP_MIX
    else if ( leqi(str,'boole-mix') ) then
       method = CC_BOOLE_MIX
    else if ( leqi(str,'mid-rule') ) then
       method = CC_MID
    else if ( leqi(str,'g-fermi') ) then
       method = CC_G_NF_0kT
       do i = G_NF_MIN_kT , G_NF_MAX_kT
          write(tmp,'(a,i0,a)') 'g-fermi(',i,')'
          if ( leqi(str,tmp) ) then
             method = CC_G_NF_0kT + i
             return
          end if
       end do
    else
       call die('Unknown method for the contour: '//trim(str))
    end if
  end function method_str

  function method_ts_c_io(c) result(ret)
    use fdf, only : leqi
    type(ts_c_io), intent(in) :: c
    integer :: ret
    ret = method(c%method)
  end function method_ts_c_io

end module m_ts_cctype
