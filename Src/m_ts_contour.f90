!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_contour
!
! Routines that are used to setup the contour for integration of the GFs
! 
! A couple of routines are included in this module
!   1) setup_contour
!   2) io_contour
!   3) print_contour
!   9) sort_contour

! Furthermore we have a couple of local routines 
! used for generating the contour points
!   4) mod_HansSkriver
!   5) sommerfeld
!   6) Gauss_Fermi_2kT_plus_line
!   7) Gauss_Fermi_0kT_plus_line
!   8) transport
!   9) phonon

! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_cctype
  use precision, only : dp

  implicit none

  ! This module will also contain all the contour variables
  integer, save :: NEn_eq
  integer, save :: NEn  ! Number of energy points in the contour path
  integer, save :: PNEn ! Number of energy points in the contour path (divisible by Nodes)

  ! Contour path
  type(ts_ccontour), dimension(:), pointer, save :: contour

  ! The contour specific variables
  real(dp), save, public :: CCEmin, GFEta, kT, EfL, EfR

  ! We need to retain the information about the contour here.
  ! It provides an easier overview as there are quite a few constants governing the
  ! methods.
  integer, save, public  :: C_Eq_Circle
  integer, save, public  :: C_Eq_Circle_N
  real(dp), save, public :: C_Eq_split
  integer, save, public  :: C_Eq_Line_bottom
  integer, save, public  :: C_Eq_Line_bottom_N
  real(dp), save, public :: C_Eq_Line_split
  integer, save, public  :: C_Eq_Line_tail
  integer, save, public  :: C_Eq_Line_tail_N
  integer, save, public  :: C_Eq_Pole_N
  integer, save, public  :: C_nEq_tail
  integer, save, public  :: C_nEq_tail_N
  real(dp), save, public :: C_nEq_split
  integer, save, public  :: C_nEq_mid
  integer, save, public  :: C_nEq_mid_N
  real(dp), save, public :: C_transport_Emin, C_transport_Emax
  integer, save, public  :: C_Transport_N

  ! The contours for the equilibrium density are attributed a fruitful discussion with
  ! Hans Skriver. Previously the routine names reflected his contribution.
  ! However, for programming clarity we have employed a different naming scheme.
  ! In order to retain the contributions it is encouraged to keep this sentiment for his
  ! contributions.

  ! Furthermore the non-equilibrium density integration are attributed discussions with
  ! Antti-Pekka Jauho. 

  ! For further attributions see the original paper by Brandbyge, et. al, 2002: DOI: 10.1103/PhysRevB.65.165401

  public :: NEn, PNEn, contour
  public :: setup_contour, io_contour, print_contour
  public :: sort_contour
  public :: contour_Eq
  public :: contour_EqL, contour_EqR, contour_nEq
  public :: contour_Transport
  public :: read_contour_options

  private

contains

  subroutine read_contour_options(kT_in, IsVolt, VoltL,VoltR)

    use units,    only : Pi, eV
    use parallel, only : IONode, Nodes, operator(.parcount.)
    use fdf

    logical, intent(in) :: IsVolt
    real(dp), intent(in) :: kT_in, VoltL, VoltR ! in Ry
    character(len=200) :: chars
    character(len=200), parameter :: OPT_CHAR = '(''ts_options: '',a,t53,''=    '',a)'
    character(len=200), parameter :: OPT_FLOAT = '(''ts_options: '',a,t53,''='',f10.4)'
    character(len=200), parameter :: OPT_INT = '(''ts_options: '',a,t53,''='',i5)'
    character(len=200), parameter :: OPT_FLOAT_UNIT = '(''ts_options: '',a,t53,''='',f10.4,tr1,a)'
    integer :: i, j, k
    real(dp) :: r
    logical :: correct

    ! Transfer options for the fermi-levels
    EfL = VoltL
    EfR = VoltR

    ! Copy over temperature
    kT = kT_in
    
    ! Read in the predefined before we can do anything else
    call fdf_deprecated('TS.ComplexContourEmin','TS.Contour.Eq.Emin')
    CCEMin = fdf_get('TS.ComplexContourEmin',-3._dp,'Ry')
    CCEMin = fdf_get('TS.Contour.Eq.Emin',CCEMin,'Ry')
    call fdf_deprecated('TS.biasContour.Eta','TS.Contour.nEq.Eta')
    GFEta = fdf_get('TS.biasContour.Eta',0.000001_dp,'Ry')
    GFEta = fdf_get('TS.Contour.nEq.Eta',GFEta,'Ry')
    if ( GFEta <= 0.d0) call die('ERROR: GFeta <= 0.0, we do not allow for &
         &using the advanced Greens function, please correct.')

    ! First we will do the defaults
    C_Eq_Circle = CC_TYPE_G_LEGENDRE
    C_Eq_Circle_N = 24
    C_Eq_split = -10._dp * kT
    ! Later down we will correct the method for the Gauss-Fermi quadrature placement
    C_Eq_Line_bottom = CC_TYPE_G_NF_0kT
    C_Eq_Line_bottom_N = 6
    C_Eq_Pole_N = 6
    ! The default bottom contour takes into account the infinite integral
    C_Eq_Line_split = 0._dp
    ! in this case the G_NF_10kT shows the same specification
    C_Eq_Line_tail = 0 ! we check against this method to determine whether the input is legit
    C_Eq_Line_tail_N = 0
    ! The non-equilibrium contour
    C_nEq_tail = CC_TYPE_G_NF_0kT - 5
    C_nEq_tail_N = 6
    C_nEq_mid = CC_TYPE_SIMP_EXT
    C_nEq_mid_N = 6

    ! ******* default setup finished ********
    ! We have now setup the default parameters for the contour method
    ! Lets read in what the user requests.

    ! First the equilibrium contour
    call fdf_deprecated('TS.ComplexContour.NPoles','TS.Contour.Eq.Pole.N')
    C_Eq_Pole_N = fdf_get('TS.ComplexContour.NPoles',C_Eq_Pole_N)
    C_Eq_Pole_N = fdf_get('TS.Contour.Eq.Pole.N',C_Eq_Pole_N)
    call fdf_deprecated('TS.ComplexContour.NCircle','TS.Contour.Eq.Circle.N')
    C_Eq_Circle_N = fdf_get('TS.ComplexContour.NCircle',C_Eq_Circle_N)
    C_Eq_Circle_N = fdf_get('TS.Contour.Eq.Circle.N',C_Eq_Circle_N)
    call fdf_deprecated('TS.ComplexContour.NLine','TS.Contour.Eq.Line.N')
    C_Eq_Line_bottom_N = fdf_get('TS.ComplexContour.NLine',C_Eq_Line_bottom_N)
    C_Eq_Line_bottom_N = fdf_get('TS.Contour.Eq.Line.N',C_Eq_Line_bottom_N)
    C_Eq_Line_bottom_N = fdf_get('TS.Contour.Eq.Line.Bottom.N',C_Eq_Line_bottom_N)
    C_Eq_Line_tail_N = fdf_get('TS.Contour.Eq.Line.Tail.N',0)

    ! The above was the standard read-in
    ! Now we read in the extra information that could be available

    ! Read in the "correct" settings for the integration
    chars = fdf_get('TS.Contour.Eq.Circle.Method','G-Legendre')
    if ( leqi(chars,'g-legendre') ) then
       C_Eq_Circle = CC_TYPE_G_LEGENDRE
    else if ( leqi(chars,'g-chebyshev-open') ) then
       C_Eq_Circle = CC_TYPE_G_CHEBYSHEV_O
    else if ( leqi(chars,'g-tschebyscheff') ) then
       C_Eq_Circle = CC_TYPE_G_TSCHEBYSHEFF
    else
       call die('Unrecognized eq. circle integration &
            &scheme: '//trim(chars))
    end if

    i = fdf_get('TS.Contour.Eq.Split.Circle-Line',-10)
    ! The splitting point between the circle and line
    C_Eq_split = i * kT

    ! Now we need to read in the method used for the bottom and tail part
    ! of the equilibrium contour
    chars = fdf_get('TS.Contour.Eq.Line.Method','g-fermi')
    i = index(chars,'+')
    if ( i == 1 ) call die('You must not start the specifier for the &
         &integration scheme with +: '//trim(chars))

    ! This will read in the tail-specification
    if ( i > 0 ) then
       i = i + 1
       ! We first figure out the tail-specification
       if ( leqi(chars(i:),'g-fermi') ) then
          ! We need to figure out how to cut it...
          ! We know that the user have specified a bottom contour method
          ! We will overwrite the user specified number seperation
          ! and determine an "optimal" splitting of the integration points
          call init_GaussFermi_Bottom_Tail(C_Eq_split, 0._dp, kT, &
               C_Eq_Line_bottom_N+C_Eq_Line_tail_N , &
               C_Eq_Line_tail_N, C_Eq_Line_bottom_N, C_Eq_Line_tail ) 
          ! Update the specified partition of the lines
          C_Eq_Line_split = C_Eq_Line_tail * kT
          C_Eq_Line_tail = CC_TYPE_G_NF_0kT + C_Eq_Line_tail

          ! when using the Gauss-Fermi method we do not allow the user
          ! to specify a split position, hence if the fdf is defined die...
          if ( fdf_defined('TS.Contour.Eq.Split.Line') ) then
             call die('Currently the Gauss-Fermi method while splitting is &
                  &fully algorithm implemented. Hence you cannot both specify: &
                  TS.Contour.Eq.Split.Line and TS.Contour.Eq.Line.Method <bot-method>+g-fermi')
          end if

          ! The default is already the Gauss-Fermi quadrature (DEFAULT)
       else if ( leqi(chars(i:),'g-laguerre') ) then
          C_Eq_Line_tail = CC_TYPE_G_LAGUERRE
       else
          call die('Unrecognized equilibrium tail integration &
               &scheme: '//trim(chars(i:)))
       end if
       i = i - 1
    end if
    
    ! Now we will read in the bottom line-segment
    i = i - 1
    if ( i <= 0 ) i = len_trim(chars)

    if ( leqi(chars(:i),'g-fermi') ) then
       ! Inifinite integrals assume that there are no specific tail integral
       ! The default is already the Gauss-Fermi quadrature (DEFAULT)
       C_Eq_Line_bottom = CC_TYPE_G_NF_0kT
    else if ( leqi(chars(:i),'g-legendre') ) then
       C_Eq_Line_bottom = CC_TYPE_G_LEGENDRE
    else if ( leqi(chars(:i),'g-chebyshev-open') ) then
       C_Eq_Line_bottom = CC_TYPE_G_CHEBYSHEV_O
    else if ( leqi(chars(:i),'g-tschebyscheff') ) then
       C_Eq_Line_bottom = CC_TYPE_G_TSCHEBYSHEFF
    else if ( leqi(chars(:i),'extended-simpson') ) then
       C_Eq_Line_bottom = CC_TYPE_SIMP_EXT
    else if ( leqi(chars(:i),'composite-simpson') ) then
       C_Eq_Line_bottom = CC_TYPE_SIMP_COMP
    else if ( leqi(chars(:i),'simpson-3/8') ) then
       C_Eq_Line_bottom = CC_TYPE_SIMP_38
    else if ( leqi(chars(:i),'mid-rule') ) then
       C_Eq_Line_bottom = CC_TYPE_MID
    else
       call die('Could not figure out the bottom contour type for the &
            &equilibrium line contour: '//trim(chars(:i)))
    end if

    select case ( C_Eq_Line_bottom ) 
    case ( CC_TYPE_G_NF_0kT )
       ! This will be catched first
       C_Eq_Line_bottom = CC_TYPE_G_NF_0kT + nint(C_Eq_split/kT)
    case ( CC_TYPE_G_NF_MIN : CC_TYPE_G_NF_0kT-1, &
         CC_TYPE_G_NF_0kT+1:CC_TYPE_G_NF_MAX )
       C_Eq_Line_bottom = CC_TYPE_G_NF_0kT + nint(C_Eq_split/kT)
    end select


    ! Determine the split on the line
    i = fdf_get('TS.Contour.Eq.Split.Line',nint(C_Eq_Line_split/kT))
    C_Eq_Line_split = i * kT

    ! Now we need to take care of how it is performed
    if ( nint(C_Eq_split / kT) > nint(C_Eq_Line_split / kT) ) then
       call die('You cannot split the circle contour before the Fermi-line &
            &integral.')
    end if
    ! General contour checks
    if ( nint(C_Eq_split / kT) >= -1 ) then
       call die('The Circle-Line split lies too close to the imaginary &
            &axis. We do not allow this due to possible overlap of the &
            &Fermi-poles and the contour. Please move the split (further) &
            &into the 2nd quadrant of the complex plane.')
    end if
    if ( CCEmin >= C_Eq_split ) then
       call die('You cannot split the circle contour before the contour &
            &beginning of the integral.')
    end if

    ! Correct for imperfect computational balance
    correct = fdf_get('TS.Contour.Eq.NoEmptyCycles',.true.)
    if ( correct ) then
       i = C_Eq_Circle_N + C_Eq_Line_bottom_N + C_Eq_Line_tail_N + C_Eq_Pole_N
       ! We immediately correct the number of energy-points for the contour
       if ( mod(i,Nodes) /= 0 ) then
          C_Eq_Circle_N = C_Eq_Circle_N + Nodes - mod(i,Nodes)
       end if
    end if

    if ( C_Eq_Line_tail /= 0 .and. C_Eq_Line_tail_N == 0 ) then
       call die('You have requested a specific tail integral for the &
            &equilibrium line. However, no points are specified for this &
            &contour. &
            &Please specify number of points with TS.Contour.Eq.Line.Tail.N')
    end if

    ! Write out how we do the integrals
    if ( IONode ) then
       write(*,opt_float_unit) 'Circle contour E_min',CCEmin / eV,'eV'
              if ( C_Eq_Circle_N > 0 ) then
          write(*,opt_char) 'Circle contour method',longtype2str(C_Eq_Circle)
       else
          call die('No points for the circle contour is defined. &
               &You need at least one point on the circle contour.')
       end if
       write(*,opt_int) 'Circle contour points',C_Eq_Circle_N
       write(*,opt_float_unit) 'Circle -> Line energy',C_Eq_split/kT,'kT'
       chars = 'Line'
       if ( C_Eq_Line_tail_N > 0 ) chars = 'Bottom line'
       if ( C_Eq_Line_bottom_N > 0 ) then
          write(*,opt_char) trim(chars)//' contour method',longtype2str(C_Eq_Line_bottom)
       else
          call die('No points for the line contour is defined. &
               &You need at least one point on the line contour.')
       end if
       
       write(*,opt_int) trim(chars)//' contour points',C_Eq_Line_bottom_N
       if ( C_Eq_Line_tail_N > 0 ) then
          write(*,opt_float_unit) 'Bottom line -> tail line energy',C_Eq_Line_split/kT,'kT'
          write(*,opt_char) 'Tail line contour method',longtype2str(C_Eq_Line_tail)
          write(*,opt_int) 'Tail line contour points',C_Eq_Line_tail_N
       end if
       if ( C_Eq_Pole_N > 0 ) then
          write(*,opt_int) 'Number of Fermi poles',C_Eq_Pole_N
       else
          call die('You need at least one pole for the Fermi function')
       end if

    end if
       
    ! *****************************************************************
    ! ****** Now we need to handle the non-equilibrium contour ********
    ! *****************************************************************

    call fdf_deprecated('TS.biasContour.NumPoints', &
         'TS.Contour.nEq.Tail.N/TS.Contour.nEq.Middle.N')

    ! Set default
    C_nEq_tail    = CC_TYPE_G_NF_0kT - 5
    C_nEq_mid     = CC_TYPE_G_NF_0kT - 5
    C_nEq_mid_N   = 6
    C_nEq_tail_N  = 6
    C_nEq_split = - 5._dp * kT

    ! Old deprecated read-ins

    if ( IsVolt ) then

       i = fdf_get('TS.biasContour.NumPoints',10)
       correct = fdf_get('TS.Contour.nEq.NoEmptyCycles',.true.)
       if ( correct .and. mod(i,Nodes) /= 0 ) i = i + Nodes - mod(i,Nodes)

       ! Determine the sizes of the individual lines
       ! Also determine the optimal Gauss-Fermi-line
       call init_GaussFermi_Bottom_Tail(0._dp, EfR, kT, i/2, &
            C_nEq_tail_N, C_nEq_mid_N, j ) 
       C_nEq_mid_N = i - 2 * C_nEq_tail_N
       C_nEq_split  = j * kT
       C_nEq_tail   = CC_TYPE_G_NF_0kT + j 

       ! Override what the user requested
       C_nEq_mid_N  = fdf_get('TS.Contour.nEq.Middle.N',C_nEq_mid_N)
       C_nEq_tail_N = fdf_get('TS.Contour.nEq.Tail.N'  ,C_nEq_tail_N)
       i = fdf_get('TS.Contour.nEq.Split',nint(C_nEq_split/kT))
       C_nEq_split  = i * kT
       C_nEq_tail   = CC_TYPE_G_NF_0kT + i

       ! The middle method is still SIMP_EXT
       
    end if

    ! The default line contour for the non-equilibrium method has already been set
    ! So we just read in the user specified methods

    ! First the default method
    call fdf_deprecated('TS.biasContour.method', 'TS.Contour.nEq.Method')
    chars = fdf_get('TS.biasContour.method','gaussfermi')
    if ( .not. fdf_defined('TS.Contour.nEq.Method') ) then
       if ( leqi(chars,'gaussfermi') ) then
          ! The gauss-fermi is already defaulted
       else if ( leqi(chars,'sommerfeld') ) then
          C_nEq_tail   = CC_TYPE_SOMMERFELD
          C_nEq_tail_N = 2 * C_nEq_tail_N + C_nEq_mid_N
          C_nEq_mid    = CC_TYPE_SOMMERFELD
          C_nEq_mid_N  = 0
       else
          call die('Could not determine contour method &
               &please use the TS.Contour keys')
       end if
    end if
    
    chars = fdf_get('TS.Contour.nEq.Method','g-fermi+extended-simpson')
    ! Figure out if there is a + in the string
    i = index(chars,'+')
    if ( i == 1 ) call die('Non-equilibrium contour method cannot be prefixed &
         &with + (then we can not determine the integration method).')
    ! Determine the middle segment method
    if ( i > 0 ) then
       i = i + 1
       if ( leqi(chars(i:),'extended-simpson') ) then
          ! This is already the default
          C_nEq_mid = CC_TYPE_SIMP_EXT
       else if ( leqi(chars(i:),'composite-simpson') ) then
          C_nEq_mid = CC_TYPE_SIMP_COMP
       else if ( leqi(chars(:i),'simpson-3/8') ) then
          C_nEq_mid = CC_TYPE_SIMP_38
       else if ( leqi(chars(i:),'mid-rule') ) then
          C_nEq_mid = CC_TYPE_MID
       else if ( leqi(chars(i:),'g-legendre') ) then
          C_nEq_mid = CC_TYPE_G_LEGENDRE
       else if ( leqi(chars(i:),'g-chebyshev-open') ) then
          C_nEq_mid = CC_TYPE_G_CHEBYSHEV_O
       else if ( leqi(chars(i:),'g-tschebyscheff') ) then
          C_nEq_mid = CC_TYPE_G_TSCHEBYSHEFF
       else
          call die('Unrecognized non-equilibrium integration &
               &scheme for the middle line: '//trim(chars(i:)))
       end if
       i = i - 1
    end if

    ! If no + is found we simulate its position
    i = i - 1
    if ( i <= 0 ) i = len_trim(chars)

    ! Determine the tail integration method
    if ( leqi(chars(1:i),'g-fermi') .or. &
         leqi(chars(1:i),'gaussfermi') ) then
       ! This is already the default and has been initialized (DEFAULT)
    else if ( leqi(chars(1:i),'g-laguerre') ) then
       C_nEq_tail = CC_TYPE_G_LAGUERRE
    else if ( leqi(chars(1:i),'sommerfeld') ) then
       C_nEq_tail = CC_TYPE_SOMMERFELD
       C_nEq_mid = CC_TYPE_SOMMERFELD
       C_nEq_tail_N = 2 * C_nEq_tail_N + C_nEq_mid_N
       C_nEq_mid_N = 0
    else
       call die('Unrecognized non-equilibrium integration &
            &scheme for the tails: '//trim(chars(1:i)))
    end if

    correct = fdf_get('TS.Contour.nEq.NoEmptyCycles',.true.)
    if ( correct ) then
       i = C_nEq_mid_N + 2 * C_nEq_tail_N
       if ( mod(i,Nodes) /= 0 ) then
          C_nEq_mid_N = C_nEq_mid_N + Nodes - mod(i,Nodes)
       end if
    end if

    ! Write out how we do the integrals
    if ( IONode .and. IsVolt ) then
       chars = 'Tail non-equilibrium contour'
       if ( C_nEq_mid_N == 0 ) chars = 'Non-equilibrium contour'
       if ( C_nEq_tail_N > 0 ) then
          write(*,opt_char) trim(chars)//' method',longtype2str(C_nEq_tail)
       else
          call die('No points for the tail non-equilibrium contour are defined. &
               &You need at least one point on the non-equilibrium contour.')
       end if
       write(*,opt_int) trim(chars)//' points',C_nEq_tail_N
       chars = 'Middle non-equilibrium contour '
       if ( C_nEq_mid_N > 0 ) then
          write(*,opt_float_unit) 'Tail -> middle contour energy',C_nEq_split/kT,'kT'
          write(*,opt_char) trim(chars)//' method',longtype2str(C_nEq_mid)
          write(*,opt_int) trim(chars)//' points',C_nEq_mid_N
       else
          call die('No points for the line contour is defined. &
               &You need at least one point on the line contour.')
       end if

    end if

    if ( IONode .and. IsVolt ) then
       i = C_Eq_Circle_N + C_Eq_Line_bottom_N + C_Eq_Line_tail_N + C_Eq_Pole_N
       if ( mod(2*i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Equilibrium energy contour points are not"
          write(*,*) "        divisable by the number of nodes."
          write(*,*) "        Better scalability is achived by changing:"
          write(*,*) "          - TS.Contour.Eq.Circle.N"
          write(*,*) "          - TS.Contour.Eq.Line.<Bottom|Tail>.N"
          write(*,*) "          - TS.Contour.Eq.Pole.N"

             ! Calculate optimal number of energy points
          i = 2*i
          write(*,'(t10,a,i4)') "Used equilibrium # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
               "Optimal equilibrium # of energy points: ",i, &
               "+- i*",Nodes
       end if
       i = C_nEq_mid_N + 2 * C_nEq_tail_N
       if ( mod(i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Non-equilibrium energy contour points are not"
          write(*,*) "        divisable by the number of nodes."
          write(*,*) "        Better scalability is achieved by changing:"
          write(*,*) "          - TS.Contour.nEq.Tail.N"
          write(*,*) "          - TS.Contour.nEq.Middle.N"
          
          ! Calculate optimal number of energy points
          write(*,'(t10,a,i4)') "Used non-equilibrium # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
               "Optimal non-equilibrium # of energy points: ",i, &
               "+- i*",Nodes
       end if
       i = C_Eq_Circle_N + C_Eq_Line_bottom_N + C_Eq_Line_tail_N + C_Eq_Pole_N
       i = 2*i + C_nEq_mid_N + 2 * C_nEq_tail_N
       if ( mod(i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Total energy contour points are not"
          write(*,*) "        divisable by the number of nodes."
          
          ! Calculate optimal number of energy points
          write(*,'(t10,a,i4)') "Used # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
               "Optimal # of energy points: ",i,"+- i*",Nodes
       end if
    else if ( IONode ) then
       i = C_Eq_Circle_N + C_Eq_Line_bottom_N + C_Eq_Line_tail_N + C_Eq_Pole_N
       
       !   - The equilibrium parts are the same computational cost
       !   * Solution make the equi contours divisible by Nodes
       if ( mod(i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Total number of energy points is &
               &not divisable by the number of nodes."
          write(*,*) "        There are no computational costs &
               &associated with increasing this."
          ! Calculate optimal number of energy points
          write(*,'(t10,a,i4)') "Used # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4)') "Optimal # of energy points: ",i
       end if
    end if

    ! Update the number of contour points
    NEn_eq = C_Eq_Circle_N + C_Eq_Line_bottom_N + C_Eq_Line_tail_N + C_Eq_Pole_N
    NEn = NEn_eq 
    if ( IsVolt ) then
       NEn = NEn_eq * 2
       NEn = NEn + C_nEq_mid_N + 2 * C_nEq_tail_N
    end if

    PNEn = Nodes .parcount. NEn

    nullify(contour)
    allocate(contour(NEn))

  end subroutine read_contour_options


  ! Routine for creating the contour
  subroutine setup_contour(IsVolt)
    
    use precision, only : dp
    use parallel,  only : IONode, Nodes, operator(.PARCOUNT.)
    use m_ts_aux, only : nf2
    use m_gauss_quad
    use m_gauss_fermi

! **********************
! * INPUT variables    *
! **********************
    logical,  intent(in) :: IsVolt ! Do we have a volt

! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c(:) => null()
    real(dp), allocatable :: x(:), w(:)
    real(dp) :: sE1, sE2, tmpE1, tmpE2
    logical :: switched
    integer :: Net, N
    integer :: i, j, k

    ! Setup all the different methods
    if ( .not. IsVolt ) then

       c => contour_Eq()
       ! We will only create the equilibrium contour
       call setup_contour_Eq(CC_PART_EQUI_CIRCLE, &
            C_Eq_Circle, C_Eq_Circle_N, C_Eq_split, &
            C_Eq_Line_bottom, C_Eq_Line_bottom_N ,C_Eq_Line_split, &
            C_Eq_Line_tail, C_Eq_Line_tail_N, C_Eq_Pole_N, &
            CCEmin, 0._dp, kT, c)

    else

       switched = EfL > EfR
       if ( switched ) then
          sE2 = EfL
          sE1 = EfR
       else
          sE1 = EfL
          sE2 = EfR
       end if

       c => contour_EqL()
       ! We will only create the equilibrium contour
       call setup_contour_Eq(CC_PART_L_EQUI_CIRCLE, &
            C_Eq_Circle, C_Eq_Circle_N, C_Eq_split, &
            C_Eq_Line_bottom, C_Eq_Line_bottom_N ,C_Eq_Line_split, &
            C_Eq_Line_tail, C_Eq_Line_tail_N, C_Eq_Pole_N, &
            CCEmin+EfL, EfL, kT, c)
       ! Note we put a minus here because the integral we want is the
       ! negative of the pole-sum and C+L integral
       c(:)%w = - c(:)%w

       c => contour_EqR()
       ! We will only create the equilibrium contour
       call setup_contour_Eq(CC_PART_R_EQUI_CIRCLE, &
            C_Eq_Circle, C_Eq_Circle_N, C_Eq_split, &
            C_Eq_Line_bottom, C_Eq_Line_bottom_N ,C_Eq_Line_split, &
            C_Eq_Line_tail, C_Eq_Line_tail_N, C_Eq_Pole_N, &
            CCEmin+EfR, EfR, kT, c)
       ! Note we put a minus here because the integral we want is the
       ! negative of the pole-sum and C+L integral
       c(:)%w = - c(:)%w

       ! The last contour is the non-equilibrium
       
       ! The voltage contour
       c => contour_nEq()
       N = size(c)
       allocate(x(N),w(N))
       Net = N+1-C_nEq_tail_N

       select case ( C_nEq_tail ) ! check the bias contour
       case ( CC_TYPE_SOMMERFELD ) ! 1. order
          if ( C_nEq_mid_N /= 0 ) call die('Sommerfeld &
               &requires no middle integration points')

          call Sommerfeld(EfL, EfR, kT, N, x, w)
          
       case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )
          if ( abs(C_nEq_split/kT-nint(C_nEq_split/kT)) > 1e-5 ) &
               call die('Splitting energy in the non-equilibrium &
               &contour is erroneous. Must be integer part of kT')

          call GaussFermi(nint(C_nEq_split/kT),C_nEq_tail_N,x(Net),w(Net))

       case ( CC_TYPE_G_LAGUERRE )
          
          call Gauss_Laguerre_GW(C_nEq_tail_N, x(Net),w(Net), &
               a=C_nEq_split)

          ! The Gauss-Laguerre is pure in the limit of 
          ! large energy, hence it makes sense to deal with it like that

          call die('not fully implemented')

       case default
          if ( IONode ) &
               write(*,*) 'ERROR: Contour not defined'
          call die('ERROR:  setup_contour: Contour not defined') 
       end select

       select case ( C_nEq_tail )
       case ( CC_TYPE_SOMMERFELD )
          ! Do nothing
          ! Type in this segment takes care of the full integral
       case default

          ! We need to symmetrize the tails
          j = C_nEq_tail_N 
          do i = Net , N
             x(j) = -x(i)
             w(j) =  w(i)
             j = j - 1
          end do

          ! This is already sorted, hence
          ! we know that we do it correct here
          sE1 = sE1 - C_nEq_split
          sE2 = sE2 + C_nEq_split

          call line_integral(C_nEq_mid,sE1, sE2, kT, C_nEq_mid_N, &
               x(C_nEq_tail_N+1), w(C_nEq_tail_N+1))
          
          sE1 = sE1 + C_nEq_split
          sE2 = sE2 - C_nEq_split

       end select

       ! We need to reverse the arguments here (positive bias is "standard")
       call correct_weight_sign(N,w,EfR,EfL)

       ! Create the tail contours
       do i = 1 , C_nEq_tail_N
          j = N + 1 - i
          tmpE1 = x(i) * kT + sE1
          tmpE2 = x(j) * kT + sE2
          c(i)%c = dcmplx(tmpE1,GFeta)
          c(j)%c = dcmplx(tmpE2,GFeta)
          select case ( C_nEq_tail ) 
          case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )
             c(i)%w = w(i) * kT
             c(j)%w = w(j) * kT
          case default
             c(i)%w = w(i) * kT * nf2(tmpE1,sE1,sE2,kT)
             c(j)%w = w(j) * kT * nf2(tmpE2,sE1,sE2,kT)
          end select
          c(i)%part = CC_PART_NON_EQUI
          c(i)%type = C_nEq_tail
          c(j)%part = CC_PART_NON_EQUI
          c(j)%type = C_nEq_tail

       end do

       ! Create the middle contours
       do i = C_nEq_tail_N + 1 , N - C_nEq_tail_N
          c(i)%c = dcmplx(x(i),GFeta)
          c(i)%w = w(i) * nf2(x(i),sE1,sE2,kT)
          c(i)%part = CC_PART_NON_EQUI
          c(i)%type = C_nEq_mid
       end do

       deallocate(x,w)
       
    end if

    
    ! Whether we should add any transport energy points for the
    ! transiesta
    if ( C_transport_N > 0 ) then
       
       c => contour_Transport()
       call transmission(C_transport_Emin, C_transport_Emax, &
            GFeta, C_transport_N, c)

    end if

  end subroutine setup_contour


! Routine for creating the contour
  subroutine setup_contour_Eq(PART, C_Eq_Circle, C_Eq_Circle_N, C_Eq_split, &
       C_Eq_Line_bottom, C_Eq_Line_bottom_N ,C_Eq_Line_split, &
       C_Eq_Line_tail, C_Eq_Line_tail_N, C_Eq_Pole_N, CCEmin, Ef, kT, contour)

    use units, only : Pi, eV
    use precision, only : dp
    use parallel,  only : IONode, Nodes, operator(.PARCOUNT.)
    use m_ts_aux, only : nf
    use m_gauss_fermi
    use m_gauss_quad

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: PART
    integer, intent(in) :: C_Eq_Circle, C_Eq_Circle_N
    integer, intent(in) :: C_Eq_Line_bottom, C_Eq_Line_bottom_N
    integer, intent(in) :: C_Eq_Line_tail, C_Eq_Line_tail_N, C_Eq_Pole_N
    real(dp), intent(in) :: C_Eq_split, C_Eq_Line_split
    real(dp), intent(in) :: CCEmin, Ef, kT

    type(ts_ccontour), pointer, intent(out) :: contour(:)
! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c(:) => null()
    real(dp), allocatable :: x(:), w(:)
    real(dp) :: Delta, D, gamma, alpha, R, beta
    complex(dp) :: z0, ztmp
    integer :: N_kT
    integer :: i, j, curN


    if ( size(contour) <= 0 ) then  ! if there are not any points we definetely have an error
       call die("ERROR: setup_contour_Eq: no energy points specified")
    end if
    
    i = 1
    curN = C_Eq_Pole_N
    if ( curN > size(contour) ) call die('Error in contour setup')
    ! Initialize the poles
    c => contour(1:curN)
    call nf_poles(PART+2,Ef,kT,C_Eq_Pole_N,c)

    ! Add the line integral (i.e. from above the last pole)
    Delta = C_Eq_Pole_N * 2._dp*Pi*kT
    
    ! Add the tail integral
    i = curN + 1
    curN = curN + C_Eq_Line_tail_N
    if ( curN > size(contour) ) call die('Error in contour setup')
    if ( C_Eq_Line_tail_N > 0 ) then
       c => contour(i:curN)
       
       ! Allocate space
       allocate(x(C_Eq_Line_tail_N),w(C_Eq_Line_tail_N))

       ! Determine the method
       select case ( C_Eq_Line_tail ) 
       case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )
          call GaussFermi(nint(C_Eq_Line_split/kT),C_Eq_Line_tail_N,x,w)

       case ( CC_TYPE_G_LAGUERRE )
          ! We retain it to be pure...
          call Gauss_Laguerre_GW(C_Eq_Line_tail_N,x,w,a=C_Eq_Line_split/kT)
          ! ...hence we need to multiply the weights by 
          w(:) = w(:)/(1._dp + exp(-x(:)-C_Eq_Line_split/kT))

       case default 
          call die('Nothing could be deciphered from the Eq. tail contour method')
       end select
    
       ! We move them over to this
       do i = 1 , C_Eq_Line_tail_N
          j = C_Eq_Line_tail_N + 1 - i
          c(i)%c    = dcmplx(x(j)*kT + Ef,Delta)
          c(i)%w    = - w(j) * kT
          c(i)%part = PART + 1
          c(i)%type = C_Eq_Line_tail
       end do
       
       deallocate(x,w)
    end if

    ! First we point to the bottom line
    i = curN + 1
    curN = curN + C_Eq_Line_bottom_N
    if ( curN > size(contour) ) call die('Error in contour setup')
    c => contour(i:curN)

    ! Allocate space
    allocate(x(C_Eq_Line_bottom_N),w(C_Eq_Line_bottom_N))

    ! Determine the method
    select case ( C_Eq_Line_bottom ) 
       ! We will not allow for -2 and 0 kT contours (maybe up to -5?)
    case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )
       if ( C_Eq_Line_tail_N > 0 ) call die('Wrong options')
       call GaussFermi(nint(C_Eq_split/kT),C_Eq_Line_bottom_N,x,w)

    case default
      if ( C_Eq_Line_tail_N == 0 ) &
           call die('Using a bounded integral requires a tail integral')

      call line_integral(C_Eq_Line_bottom, C_Eq_split/kT, &
           C_Eq_Line_split/kT, kT, C_Eq_Line_bottom_N,x,w)
      w(:) = w(:) * nf(x(:) * kT + Ef)
       
   end select

    ! We move them over to this
    do i = 1 , C_Eq_Line_bottom_N
       j = C_Eq_Line_bottom_N + 1 - i
       c(i)%c    = dcmplx(x(j)*kT + Ef,Delta)
       c(i)%w    = - w(j) * kT
       c(i)%part = PART + 1
       c(i)%type = C_Eq_Line_bottom
    end do

    deallocate(x,w)

    ! We create the circle contour
    D     = Ef - CCEmin
    ! The ending point for the contour (on the imaginary axis)
    gamma = - C_Eq_split
    ! The angle between CCEmin and the point)
    alpha = dATAN(Delta/(D-gamma))
    ! The radius of the circle in the complex plane
    R     = dsqrt(Delta**2 + (D - gamma)**2) / (2._dp*cos(alpha))
    ! The center of the circle on the real-axis
    z0    = dcmplx(CCEmin + R, 0._dp)
    ! The angle from where we start
    beta  = dasin(Delta / R)

    ! First we point to the bottom line
    i = curN + 1
    curN = curN + C_Eq_Circle_N
    if ( curN > size(contour) ) call die('Error in contour setup')
    c => contour(i:curN)
    
    ! Allocate space
    allocate(x(C_Eq_Circle_N),w(C_Eq_Circle_N))

    select case ( C_Eq_Circle )
    case ( CC_TYPE_G_LEGENDRE )
       call Gauss_Legendre_Rec(C_Eq_Circle_N, 0, beta, Pi, x, w)

    case ( CC_TYPE_G_CHEBYSHEV_O )
       ! We switch the integration bounds and change the sign of the weight
       ! this will give the correct path when out-put
       call Gauss_Chebyshev_Exact(C_Eq_Circle_N, x, w, a=beta, b=Pi, &
            method=1, pure=.false.)

    case ( CC_TYPE_G_CHEBYSHEV_C )
       ! We switch the integration bounds and change the sign of the weight
       ! this will give the correct path when out-put
       call Gauss_Chebyshev_Exact(C_Eq_Circle_N, x, w, a=beta, b=Pi, &
            method=0, pure=.false.)

    case ( CC_TYPE_G_TSCHEBYSHEFF )
       D = Pi - beta
       call Gauss_Tschebyscheff_Exact(C_Eq_Circle_N,x,w, &
            b=D,pure=.false.)
       x(:) = x(:) + beta

    case default
       call die('Unrecognized contour for the equilibrium circle')
    end select

    do i = 1 , C_Eq_Circle_N 
       ztmp   = R * cdexp(dcmplx(0._dp,x(i)))
       c(i)%c = z0 + ztmp
       ! Factor i, comes from Ed\theta=dE=iR e^{i\theta}
       c(i)%w = w(i) * nf((c(i)%c-Ef)/kT) * dcmplx(0._dp,1._dp)*ztmp
       c(i)%part = PART
       c(i)%type = C_Eq_Circle
    end do

    deallocate(x,w)
    
  end subroutine setup_contour_Eq

!------------------------------------------------------------
!
!     1. order Sommerfeld expansion using 2kT increments
!      
! We assume that the function to be integrated varies slowly on the
! kT-scale         
  subroutine Sommerfeld(E1,E2, kT, NEn,x,w)

    use precision, only : dp
    use parallel, only : IONode
    use units, only : Pi

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E1, E2        ! energy parameters 
    real(dp), intent(in) :: kT            ! temperature in Ry
    integer,  intent(in) :: NEn

! ***********************
! * OUTPUT variables    *
! ***********************
    real(dp), intent(out) :: x(NEn), w(NEn)

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: EE1, EE2
    real(dp) :: delta, etaSF, rtmp
    integer :: i

    if ( NEn <= 4 ) then
       if(IONode) &
            write(*,*) &
            'ERROR: No. points=',NEn,' not valid for Sommerfeld'
       call die('ERROR: Contour: no. points not OK for Sommerfeld') 
    end if
    
    EE2 = max(E1,E2)
    EE1 = min(E1,E2)
       
    delta = (EE2 - EE1)/real(NEn-3,dp) 
       
    etaSF = kT
    
    rtmp = kT*kT*Pi*Pi/(12._dp*etaSF)
    
    ! Assign the contour points and weights
    x(1) = EE1 - etaSF
    w(1) = 0.25_dp * delta + rtmp
    x(2) = EE1 + etaSF
    w(2) = 0.25_dp * delta - rtmp
    
    x(NEn-1) = EE2 - etaSF
    w(NEn-1) = 0.25_dp * delta - rtmp
    
    x(NEn)   = EE2 + etaSF
    w(NEn)   = 0.25_dp * delta + rtmp
    
    do i = 3 , NEn - 2
       x(i) = delta * (i-2) + EE1
       w(i) = delta
    end do
     
  end subroutine Sommerfeld

  subroutine line_integral(TYPE,EE1,EE2,kT,NEn,x,w)
    
    use precision, only : dp
    use units, only : Pi
    use m_ts_cctype
    use m_gauss_quad ! Just all, many routines can be used

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: TYPE
    real(dp), intent(in) :: EE1,EE2    ! energy parameters
    real(dp), intent(in) :: kT       ! temperature in Ry
    integer,  intent(in) :: NEn      ! No. contour points

! ***********************
! * OUTPUT variables    *
! ***********************
    real(dp), intent(out) :: x(NEn), w(NEn)

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: E1,E2 ! Sorted energies
    integer :: i, Ni
    real(dp) :: delta
    
    ! sort the energies
    E1 = min(EE1,EE2)
    E2 = max(EE1,EE2)
    
    select case ( TYPE )
    case ( CC_TYPE_SIMP_EXT )
       if ( NEn < 8 ) call die('Cannot do extended Simpson on these &
            &points. Please use more points or another method.')
       ! in simpson we count: i = 0, ..., N
       Ni = NEn - 1
       
       delta = (E2 - E1)/real(Ni,dp)
       do i = 1 , NEn
          x(i) = E1 + delta * (i-1)
          w(i) = delta
       end do

       ! extended simpsons rule
       w(1    ) = w(1    )*17.d0/48.d0
       w(2    ) = w(2    )*59.d0/48.d0
       w(3    ) = w(3    )*43.d0/48.d0
       w(4    ) = w(4    )*49.d0/48.d0
       w(NEn-3) = w(NEn-3)*49.d0/48.d0
       w(NEn-2) = w(NEn-2)*43.d0/48.d0
       w(NEn-1) = w(NEn-1)*59.d0/48.d0
       w(NEn  ) = w(NEn  )*17.d0/48.d0

    case ( CC_TYPE_SIMP_COMP )

       ! in simpson we count: i = 0, ..., N
       Ni = NEn - 1
       if ( mod(Ni,2) /= 0 ) call die('Composite Simpson rule &
            &is only allowed for uneven N, please increase points.')
       
       delta = (E2 - E1)/real(Ni,dp)
       do i = 1 , NEn
          x(i) = E1 + delta * (i-1)
          w(i) = delta / 3._dp
       end do
       
       ! Correct the weights for the composite simpson rule
       do i = 2 , NEn - 1, 2
          w(i) = w(i) * 4._dp
       end do
       do i = 3 , NEn - 2, 2
          w(i) = w(i) * 2._dp
       end do

    case ( CC_TYPE_SIMP_38 )

       ! in simpson we count: i = 0, ..., N
       Ni = NEn - 1
       if ( mod(Ni,2) /= 0 ) call die('Simpson 3/8 rule &
            &is only allowed for uneven N, please increase points.')

       delta = (E2 - E1)/real(Ni,dp)
       do i = 1 , NEn
          x(i) = E1 + delta * (i-1)
          w(i) = delta * 3._dp / 8._dp
       end do

       do i = 2 , NEn - 1, 3
          w(i)   = w(i)   * 3._dp
          w(i+1) = w(i+1) * 3._dp
          w(i+2) = w(i+2) * 2._dp
       end do
       
    case ( CC_TYPE_MID )

       ! set boundaries for gaussian quadrature
       delta = (E2 - E1)/real(NEn,dp)
       do i = 1 , NEn
          ! move into the middle of the current segment
          x(i) = E1 + delta * ( real(i,dp)-.5_dp )
          w(i) = delta
       end do

    case ( CC_TYPE_LEFT )

       ! set boundaries for gaussian quadrature
       delta = (E2 - E1)/real(NEn,dp)
       do i = 1 , NEn
          ! move into the middle of the current segment
          x(i) = E1 + delta * real(i-1,dp)
          w(i) = delta
       end do

    case ( CC_TYPE_RIGHT )

       ! set boundaries for gaussian quadrature
       delta = (E2 - E1)/real(NEn,dp)
       do i = 1 , NEn
          ! move into the middle of the current segment
          x(i) = E1 + delta * real(i,dp)
          w(i) = delta
       end do

    case ( CC_TYPE_G_LEGENDRE ) 

       call Gauss_Legendre_Rec(NEn, 0, E1, E2, x, w)

    case ( CC_TYPE_G_CHEBYSHEV_O ) 

       call Gauss_Chebyshev_Exact(NEn, x, w, &
            a=E1, b=E2, method=1, pure=.false.)

    case ( CC_TYPE_G_CHEBYSHEV_C ) 

       call Gauss_Chebyshev_Exact(NEn, x, w, &
            a=E1, b=E2, method=0, pure=.false.)

    case ( CC_TYPE_G_TSCHEBYSHEFF ) 

       delta = E2 - E1
       call Gauss_Tschebyscheff_Exact(NEn, x, w, &
            b=delta, pure=.false.)
       x(:) = x(:) + E1

       call die('Not fully implemented')

    case default
       call die('Could not determine the non-equilibrium line contour')
    end select

  end subroutine line_integral
  

  ! The residuals of the fermi-function at a real-energy
  subroutine nf_poles(PART,E, kT, Npol,c)

    use precision, only : dp
    use units, only : Pi

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: PART ! part of the contour
    real(dp), intent(in) :: E    ! at energy
    real(dp), intent(in) :: kT   ! temperature in Ry
    integer,  intent(in) :: Npol

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), target, intent(out) :: c(Npol)

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: i

    ! Residuals
    do i = 1 , Npol 
       c(i)%c    = dcmplx(E    , Pi*kT*(2._dp*(i-1)+1._dp))
       c(i)%w    = dcmplx(0._dp, Pi*kT* 2._dp) 
       c(i)%part = PART
       c(i)%type = CC_TYPE_RES
    end do
    
  end subroutine nf_poles


  
  ! Determine the number of tail-points in the Gauss-Fermi curve
  subroutine init_GaussFermi_Bottom_Tail(E1, E2, &
       kT, NEn, Ntail, Nmid, GF_N_kT)
    
    use units, only : Pi
    use precision, only : dp
    use parallel, only : IONode
    use m_gauss_fermi

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E1, E2  ! energy parameters 
    real(dp), intent(in) :: kT      ! temperature in Ry
    integer,  intent(in) :: NEn

! ***********************
! * OUTPUT variables    *
! ***********************
    integer,  intent(out) :: Ntail, Nmid, GF_N_kT

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: deltaN, VkT

    VkT = dabs(E2-E1) / kT

    ! Determine the optimal Gauss-Fermi curve
    ! If the bias is larger than 10 * kT then we can select from [-5kT;inf]
    ! almost without-overlap
    if ( VkT > 5._dp ) then
       ! we default to select -5kT (the fermi function is almost zero there)
       GF_N_kT = -5
    else
       ! We need to figure out were it makes sense to cut them
       ! Calculate the overlap of the fermi functions
       GF_N_kT = nint((5._dp-VkT)*.5_dp + .5_dp)
    end if

    ! Select an allowed range (this should not be a problem)
    GF_N_kT = min(GF_N_kT,G_NF_MAX_kT)
    GF_N_kT = max(GF_N_kT,G_NF_MIN_kT)
       
    ! Now we need to determine the number of points in each segment
    ! We do this by assuming that the weight in the tails equals that of the overlap
    ! region of the Fermi function.
    ! This seems like a reasonable assumption.

    ! Calculate the parts in the tail
    deltaN = VkT / real(NEn,dp)

    if ( GF_N_kT > 0 ) then
       ! When extending the integral into the positive region
       ! it will be better to allow the full range to be divided
       deltaN = (VkT + GF_N_kT) / real(NEn,dp)
       Ntail = nint(GF_N_kT/deltaN)
    else
       Ntail = nint(real(abs(GF_N_kT),dp)/deltaN)
    end if
    
    ! correct for spill-out
    Ntail = min(Ntail,G_NF_MAX_N)
    Ntail = max(Ntail,G_NF_MIN_N)

    ! Calculate the number of points in the middle segment
    Nmid  = NEn - Ntail

  end subroutine init_GaussFermi_Bottom_Tail
  

! ##################################################################
! ##      Generate (close to) real axis energy contour            ##
! ##           Transmission and DOS calculation                   ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
  subroutine transmission(E1,E2,GFeta,NEn,contour)

    use precision, only : dp

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E1,E2    ! energy parameters 
    real(dp), intent(in) :: GFeta    ! state broadening in Ry
    integer,  intent(in) :: NEn      ! No. contour points

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), intent(out) :: contour(NEn)   ! points for contour

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: delta
    integer :: ic

    ! Calculate the weight of each contour point
    delta = (E2-E1)/(1.d0*max(NEn-1,1))

    do ic = 1 , NEn
       contour(ic)%c = dcmplx(E1+(ic-1)*delta, GFeta)
       contour(ic)%w = dcmplx(delta          , 0d0)
       contour(ic)%part = CC_PART_TRANSPORT
       ! ???
       ! The part is a transport point, so perhaps dubious to 
       ! also include a type?
       ! ???
       ! Right now we use it for printing purposes...
       contour(ic)%type = CC_TYPE_TRANSPORT
    end do

  end subroutine transmission


! ##################################################################
! ##      Generate a phonon energy contour to be able to          ##
! ##                 handle phonon transport                      ##
! ##           Transmission and DOS calculation                   ##
! ##                            By                                ##
! ##        Nick Papior Andersen, nickpapior@gmail.com            ##
! ##################################################################
  subroutine phonon(E1,E2,GFeta,NEn,contour)

    use precision, only : dp

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E1,E2    ! energy parameters 
    real(dp), intent(in) :: GFeta    ! state broadening in Ry
    integer,  intent(in) :: NEn      ! No. contour points

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), intent(out) :: contour(NEn)   ! points for contour

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: delta
    integer :: ic

    ! Calculate the weight of each contour point
    delta = (E2-E1)/(1.d0*max(NEn-1,1))

    if ( E1 <= 0._dp .or. E2 <= 0._dp ) then
       call die('Energy range for phonon transport must not &
            &coincide with: E <= 0')
    end if

    do ic = 1 , NEn
       ! The energy contour in Phonon space is:
       ! (\omega + i \eta)**2
       contour(ic)%c = dcmplx(E1+(ic-1)*delta, GFeta)**2
       contour(ic)%w = dcmplx(delta          , 0d0)
       contour(ic)%part = CC_PART_TRANSPORT
       ! ???
       ! The part is a transport point, so perhaps dubious to 
       ! also include a type?
       ! ???
       ! Right now we use it for printing purposes...
       contour(ic)%type = CC_TYPE_TRANS_PHONON
    end do

  end subroutine phonon


  ! In order to retain the best numerical accuracy we introduce
  ! an ordering of the energy contour by weights.
  ! It is the only reasonable thing to access as the functional is
  ! non-deterministic.
  ! An example of the importance of this sorting:
  ! Take +200 voltage energy points. 
  ! The weights will be something like this:
  ! 1e-5,...,1e-4,...,1e-3,...,1e-4,...,1e-5
  ! Which means that the summation up till the half works great.
  ! But when we reach the last weight we could be in the situation where:
  !  1._dp + 1e-12_dp which will limit the accuracy obtained for that energy point.

  ! In order to circumvent this we simply sort by weights.
  subroutine sort_contour(NC,c)
    integer, intent(in) :: NC
    type(ts_ccontour), intent(inout) :: c(NC)
    
    ! Local variables
    type(ts_ccontour) :: ctmp
    integer :: i,j
    ! As we will only do this once we dont need a fancy sorting
    ! algorithm...
    do i = 1 , NC - 1
       ctmp = c(i)
       do j = i+1, NC
          if ( real(ctmp%w*conjg(ctmp%w)) > &
               real(c(j)%w*conjg(c(j)%w)) ) then
             c(i) = c(j)
             c(j) = ctmp
             ctmp = c(i)
          end if
       end do
    end do

  end subroutine sort_contour

  ! Routine for "pretty" printing the contour points out in the out file
  subroutine print_contour()
    use parallel, only : IONode
    use units,    only : eV

! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c
    character(len=CC_TYPE_LEN) :: ctype
    character(len=3) :: f_type
    integer :: i, part, type
    
    ! Initialize variables
    part = -1
    type = -1
    nullify(c)

    if ( IONode ) then
       write(*,'(a)') "transiesta: contour integration path:"
       write(f_type,'(a,i2)') 'a',CC_TYPE_LEN
       write(*,'(1x,'//trim(f_type)//',''   '',2(tr1,a12),2(tr1,a14))') &
            "Type  ","Re(c)[eV]","Im(c)[eV]","Re(weight)","Im(weight)"
       do i = 1 , NEn
          ! loop !
          c => contour(i)
          if ( part /= c%part ) then
             write(*,'(1x,a)') part2str(c)
             part = c%part
          end if
          if ( type /= c%type ) then
             ctype = type2str(c)
             type = c%type
          end if
          
         ! Write out the contour information:
          write(*,'(1x,'//trim(f_type)//','' : '',tr1,2(f12.5,tr1),2(f14.9,tr1))') &
               ctype,c%c/eV,c%w
       end do
       write(*,*) ! New line
    end if
  end subroutine print_contour


! Write out the contour to a contour file
  subroutine io_contour(slabel)
    use parallel, only : IONode
    use units, only : eV
    character(len=*), intent(in) :: slabel

! *********************
! * LOCAL variables   *
! *********************
    character(len=len_trim(slabel)+9) :: fname
    integer :: i, unit, part

    interface
       function paste(s1,s2)
         character(LEN=*), intent(in) :: s1,s2
         character(LEN=200) :: paste
       end function paste
    end interface

    if ( IONode ) then
       fname = trim(paste(slabel,'.TSCC'))
       call io_assign( unit )
       open( unit, file=fname, status='unknown' )
       write(unit,'(a)') "# Complex contour path"
       write(unit,'(4(a12,tr1))') "# Re(c)[eV] ","Im(c)[eV]","Re(w)","Im(w)"
       part = contour(1)%part
       do i = 1 , NEn
          if ( part /= contour(i)%part ) then
             write(unit,*)
             part = contour(i)%part
          end if
          write(unit,'(4(f12.6,tr1))') contour(i)%c/eV,contour(i)%w
       end do
       call io_close( unit )
    end if

  end subroutine io_contour

  subroutine correct_weight_sign(N,w,E1,E2)
    integer, intent(in) :: N
    real(dp), intent(inout) :: w(N)
    real(dp), intent(in) :: E1,E2
    if ( E1 > E2 ) then
       w(:) = -w(:)
    end if
  end subroutine correct_weight_sign

  function contour_Eq() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour_EqL()
  end function contour_Eq

  function contour_EqL() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(1:NEn_eq)
  end function contour_EqL

  function contour_EqR() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(NEn_eq+1:NEn_eq*2)
  end function contour_EqR

  function contour_nEq() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(NEn_eq*2+1:NEn-C_Transport_N)
  end function contour_nEq

  function contour_Transport() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(NEn-C_Transport_N:NEn)
  end function contour_Transport

end module m_ts_contour
