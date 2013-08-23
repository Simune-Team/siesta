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

module m_ts_io_contour
!
! Routines that are used to read in and print out the contour for integration of the GFs
! 
  use m_ts_cctype
  use precision, only : dp

  implicit none

  ! The character length of the contour text...
  integer, parameter :: c_ab_N = 200
  ! We have a type for reading in the contour and defines the contour type segments
  type :: ts_io_c
     sequence
     ! The integral bounds
     real(dp) :: a, b
     ! Number of points in the integral
     integer :: N
     ! the method employed
     integer :: type
     ! optional arguments (currently only one optional argument can be specified)
     ! In case of Gauss-Fermi this is the E_F of the contour segment...
     real(dp) :: opt
     ! optional arguments (currently only one optional argument can be specified)
     real(dp) :: G_F_off
     integer  :: G_F_lower, G_F_upper
     character(len=c_ab_N) :: c_a, c_b
  end type ts_io_c

contains

  subroutine ts_read_contour_options(cEq,cnEq, kT, IsVolt, V_half)

    use parallel, only : IONode, Nodes, operator(.parcount.)
    use fdf

    type(ts_io_c), intent(inout), allocatable :: cEq(:), cnEq(:)
    real(dp), intent(in) :: kT, V_half ! in Ry
    logical :: IsVolt
    character(len=200) :: chars
    integer :: i,j
    logical :: correct

    ! we setup the default contour types...
    if ( allocated(cEq) ) deallocate(cEq)
    allocate(cEq(2))
    if ( allocated(cnEq) ) deallocate(cnEq)
    allocate(cnEq(3))

    ! *****************************************************************
    ! ****** Now we need to handle the equilibrium contour     ********
    ! *****************************************************************
    
    ! This is the circle contour...
    cEq(1)%type = CC_TYPE_G_LEGENDRE
    cEq(1)%N    = 24
    cEq(1)%a    =  -3._dp
    cEq(1)%b    = -10._dp * kT
    cEq(1)%c_a  = '-3 Ry'
    cEq(1)%c_b  = '-10 kT'

    ! This is the line contour...
    cEq(2)%type = CC_TYPE_G_NF_0kT - 10
    cEq(2)%N    = 6
    cEq(2)%a    = cEq(1)%b
    cEq(2)%b    = huge(0._dp)
    cEq(2)%c_a  = '-10 kT'
    cEq(2)%c_b  = 'inf'
    cEq(2)%G_F_off = 0._dp
    cEq(2)%G_F_lower = -10
    cEq(2)%G_F_upper = huge(1)

    ! Lets read in what the user requests.
    call fdf_deprecated('TS.ComplexContourEmin','TS.Contour.Eq.Emin')
    cEq(1)%a = fdf_get('TS.ComplexContourEmin',cEq(1)%a,'Ry')
    cEq(1)%a = fdf_get('TS.Contour.Eq.Emin',cEq(1)%a,'Ry')
    ! ensure that it is written to the string (currently this will always be
    ! in Ry)
    write(cEq(1)%c_a,'(g20.13,tr1,a)') cEq(1)%a,'Ry'
    call correct_E_string(cEq(1)%c_a)

    call fdf_deprecated('TS.ComplexContour.NCircle','TS.Contour.Eq.Circle.N')
    cEq(1)%N = fdf_get('TS.ComplexContour.NCircle',cEq(1)%N)
    cEq(1)%N = fdf_get('TS.Contour.Eq.Circle.N',cEq(1)%N)
    call fdf_deprecated('TS.ComplexContour.NLine','TS.Contour.Eq.Line.N')
    cEq(2)%N = fdf_get('TS.ComplexContour.NLine',cEq(2)%N)
    cEq(2)%N = fdf_get('TS.Contour.Eq.Line.N',cEq(2)%N)

    ! The above was the standard read-in
    ! Now we read in the extra information that could be available
       
    ! *****************************************************************
    ! ****** Now we need to handle the non-equilibrium contour ********
    ! *****************************************************************

    if ( IsVolt ) then

    ! read in number of points for the non-equilibrium contour.
    call fdf_deprecated('TS.biasContour.NumPoints', 'TS.Contour.nEq.N')
    i = fdf_get('TS.biasContour.NumPoints',10)
    i = fdf_get('TS,Contour.nEq.N',i)

    if ( i <= 0 .and. abs(V_half) > 0._dp ) then
       call die('You seem to have requested zero bias contour &
            &points in a bias calculation. Please correct input.')
    end if

    ! Determine the sizes of the individual lines
    ! Also determine the optimal Gauss-Fermi-line
    call init_GaussFermi(0._dp, abs(V_half), kT, i/2, &
         cnEq(2)%N, cnEq(1)%N, cnEq(1)%G_F_lower ) 

    ! The left tail integral (N is set)
    cnEq(1)%type = CC_TYPE_G_NF_0kT + cnEq(1)%G_F_lower
    cnEq(1)%a = -huge(0._dp)
    cnEq(1)%b = -abs(V_half) - kT * cnEq(1)%G_F_lower
    cnEq(1)%c_a = '-inf'
    cnEq(1)%G_F_off = -abs(V_half)
    cnEq(1)%G_F_upper = huge(1)
    if ( cnEq(1)%G_F_lower > 0 ) then
       write(cnEq(1)%c_b,'(''- V/2 - '',i0,'' kT'')') abs(cnEq(1)%G_F_lower)
    else
       write(cnEq(1)%c_b,'(''- V/2 + '',i0,'' kT'')') abs(cnEq(1)%G_F_lower)
    end if

    ! The right tail integral
    cnEq(3)%type = CC_TYPE_G_NF_0kT + cnEq(1)%G_F_lower
    cnEq(3)%a = abs(V_half) + kT * cnEq(1)%G_F_lower
    cnEq(3)%b = huge(0._dp)
    cnEq(3)%N = cnEq(1)%N ! tails are the same...
    cnEq(3)%G_F_off = abs(V_half)
    cnEq(3)%G_F_lower = cnEq(1)%G_F_lower
    cnEq(3)%G_F_upper = huge(1)
    if ( cnEq(1)%G_F_lower < 0 ) then
       write(cnEq(3)%c_a,'(''V/2 - '',i0,'' kT'')') abs(cnEq(1)%G_F_lower)
    else
       write(cnEq(3)%c_a,'(''V/2 + '',i0,'' kT'')') abs(cnEq(1)%G_F_lower)
    end if
    cnEq(3)%c_b = 'inf'
    
    ! The middle integral (N is set)
    cnEq(2)%type = CC_TYPE_SIMP_MIX
    cnEq(2)%a = cnEq(1)%b
    cnEq(2)%b = cnEq(3)%a
    cnEq(2)%c_a = 'prev'
    cnEq(2)%c_b = 'next'

    ! ensure that the specified number of points
    ! is used...
    j = i - sum(cnEq(:)%N)
    if ( j /= 0 ) then
       cnEq(2)%N = cnEq(2)%N + j
       if ( cnEq(2)%N <= 0 ) then
          call die('Please increase number of bias points. &
               &Or input manually the entire block. &
               &We cannot determine the optimal split of the Gauss-Fermi &
               &function.')
       end if
    end if

    ! The default line contour for the non-equilibrium
    ! method has already been set
    ! So we just read in the user specified methods

    ! First the default method
    call fdf_deprecated('TS.biasContour.method', 'TS.Contour.nEq.Method')
    chars = fdf_get('TS.biasContour.method','g-fermi')
    chars = fdf_get('TS.Contour.nEq.Method',trim(chars))
    if ( leqi(chars,'gaussfermi') .or. &
         leqi(chars,'g-fermi') ) then
       ! The gauss-fermi is already defaulted
    else if ( leqi(chars,'sommerfeld') ) then
       i = sum(cnEq(:)%N)
       deallocate(cnEq)
       allocate(cnEq(1))
       cnEq(1)%type = CC_TYPE_SOMMERFELD
       cnEq(1)%N    = i
       ! this really doesn't matter, however, we set
       ! them anyway...
       cnEq(1)%a    = -huge(0._dp)
       cnEq(1)%b    =  huge(0._dp)
       cnEq(1)%c_a  = '-inf'
       cnEq(1)%c_b  =  'inf'
    else
       call die('Could not determine contour method: '// &
            trim(chars)//'. Please use the TS.Contour.nEq block.')
    end if

    end if
    
    ! *****************************************************************
    ! ****** Now we need to handle the block contours          ********
    ! *****************************************************************
    ! Remember that the block-contours has precedence for the setup
    ! of the contour.
    ! (return value 'i' is merely the size of cEq/cnEq)
    call ts_read_contour_block('TS.Contour.Eq' ,i,cEq ,kT)
    call check_Gauss_Fermi('TS.Contour.Eq',cEq ,kT)
    if ( IsVolt ) then
       call ts_read_contour_block('TS.Contour.nEq',i,cnEq,kT,V_half)
       call check_Gauss_Fermi('TS.Contour.nEq',cnEq,kT)
    end if

    ! We need to check whether the segments are correct
    ! The restrictions checked for can be seen in the die commands...
    j = 1
    correct = .true.
    do i = 1 , size(cEq) - 1
       ! We check that the segments are connected
       if ( dabs(cEq(i)%b-cEq(i+1)%a) > 1.e-6_dp ) then
          chars = 'The equilibrium contour segments are not connected. &
               &Please ensure that it is a continuous integral.'
          correct = .false.
       end if
       select case ( cEq(i)%type )
       case ( CC_TYPE_G_LEGENDRE, CC_TYPE_TANH_SINH, CC_TYPE_SIMP_MIX )
          ! fine
       case ( CC_TYPE_BOOLE_MIX, CC_TYPE_MID )
          ! fine
       case default
          chars = 'The equilibrium circle and middle segment contours &
               &are only allowed to use quadrature/Newton-Cotes methods &
               &with weight function w(x) = 1. Please restrict the input &
               &to comply with this.'
          correct = .false.
       end select
    end do
        
    ! The last segment of the equilibrium MUST be a Gauss-Fermi contour
    i = size(cEq)
    select case ( cEq(i)%type )
    case ( CC_TYPE_G_NF_MIN : CC_TYPE_G_NF_MAX )
       ! perfect
    case default
       chars = 'The tail integral of the Equilibrium contour one &
            &must specify the Gauss-Fermi integral. Currently it is &
            &the only supported quadrature.'
       correct = .false.
    end select

    ! in case we have a bias
    if ( IsVolt ) then
       ! The non-equilibrium tails MUST be Gauss-Fermi... (or a Sommerfeld)
       i = size(cnEq)
       if ( i == 1 ) then
          if ( cnEq(i)%type == CC_TYPE_SOMMERFELD ) then
             ! fine
          else
             chars = 'The non-equilibrium contour must be a Sommerfeld &
                  &contour in case of only one segment. &
                  &Please add segments or use the Sommerfeld.'
             correct = .false.
             j = 2
          end if
       else
          ! We check that the segments are connected
          if ( dabs(cnEq(1)%b-cnEq(2)%a) > 1.e-6_dp ) then
             chars = 'The non-equilibrium contour segments are not connected. &
                  &Please ensure that it is a continuous integral.'
             correct = .false.
             j = 2
          end if

          select case ( cnEq(1)%type )
          case ( CC_TYPE_G_NF_MIN : CC_TYPE_G_NF_MAX )
             ! perfect
          case default
             chars = 'The tail integral of the non-equilibrium contour one &
                  &must specify the Gauss-Fermi integral. Currently it is &
                  &the only supported quadrature.'
             correct = .false.
             j = 2
          end select
          select case ( cnEq(i)%type )
          case ( CC_TYPE_G_NF_MIN : CC_TYPE_G_NF_MAX )
             ! perfect
          case default
             chars = 'The tail integral of the non-equilibrium contour one &
                  &must specify the Gauss-Fermi integral. Currently it is &
                  &the only supported quadrature.'
             correct = .false.
             j = 2
          end select

          do i = 2 , size(cnEq) - 1
             ! We check that the segments are connected
             if ( dabs(cnEq(i)%b-cnEq(i+1)%a) > 1.e-6_dp ) then
                chars = 'The non-equilibrium contour segments are not connected. &
                     &Please ensure that it is a continuous integral.'
                correct = .false.
                j = 2
             end if

             select case ( cnEq(i)%type )
             case ( CC_TYPE_G_LEGENDRE, CC_TYPE_TANH_SINH, CC_TYPE_SIMP_MIX )
                ! fine
             case ( CC_TYPE_BOOLE_MIX, CC_TYPE_MID )
                ! fine
             case default
                chars = 'The non-equilibrium middle segment contours &
                     &are only allowed to use quadrature/Newton-Cotes methods &
                     &with weight function w(x) = 1. Please restrict the input &
                     &to comply with this.'
                correct = .false.
                j = 2
             end select
          end do
       end if
       
    end if

    if ( .not. correct ) then
       if ( j == 1 ) then
          call ts_print_contour_block('TS.Contour.Eq',cEq)
       else
          call ts_print_contour_block('TS.Contour.nEq',cnEq)
       end if
       call die(trim(chars))
    end if

 
  contains

    subroutine check_Gauss_Fermi(bName,c,kT)
      character(len=*), intent(in) :: bName
      type(ts_io_c), intent(in) :: c(:)
      real(dp), intent(in) :: kT
      logical :: error
      integer :: i, ikT

      error = .false.
      do i = 1 , size(c)
         select case ( c(i)%type )
         case ( G_NF_MIN_kT : G_NF_MAX_kT )
            ! We need to ensure that the end-point integral
            ! is well defined...
            if ( abs(c(i)%a) < abs(c(i)%b) ) then
               ! We have an integral in the positive end
               if ( abs(c(i)%b) < huge(0._dp) ) then
                  ikT = nint(abs(c(i)%b) / kT)
                  if ( ikT /= c(i)%G_F_upper ) call die('An inconsistency &
                       &please contact the developer')
               end if
            else
               if ( abs(c(i)%a) < huge(0._dp) ) then
                  ikT = nint(abs(c(i)%a) / kT)
                  if ( ikT /= c(i)%G_F_upper ) call die('An inconsistency &
                       &please contact the developer')
               end if
            end if
            select case ( ikT )
            case ( 17:20,22,24,26,28,30 )
               ! do nothing, these are fine
            case default
               call die('The Gauss-Fermi contour only allow &
                    &17, 18, 19, 20, 22, 24, 26, 28 or 30 &
                    &as input, for infinite integral write inf.')
            end select
         end select

         ! we do a check of the contours if they are error-prone, we will write
         ! the interpretation out and quit.
         ! the most important thing to check is the ascending order:

         ! Check for ascending order
         error = error .or. c(i)%a >= c(i)%b
         if ( i > 1 ) then
            error = error .or. ( c(i-1)%b > c(i)%a )
            if ( leqi(c(i-1)%c_b,'next') ) then
               error = error .or. ( c(i-1)%a >= c(i-1)%b )
            end if
         end if
         
      end do

      if ( error ) then
         ! write out the block as interpreted....
         call ts_print_contour_block(bName,c)

         if ( IONode ) then
            write(*,'(a)') 'Values in the bounds lines (same order as block)'
            do i = 1 , size(c)
               ! print out the values
               write(*,'(t5,g20.13,a,g20.13)') c(i)%a,'   to   ',c(i)%b
            end do
         end if
         
         call die('The integrals must be in ascending order. Please &
              &correct block: '//trim(bName))
      end if

    end subroutine check_Gauss_Fermi

    ! Determine the number of tail-points in the Gauss-Fermi curve
    subroutine init_GaussFermi(E1, E2, kT, &
         N, N_mid, N_tail, GF_N_kT)

      use units, only : Pi
      use precision, only : dp
      use parallel, only : IONode
      use m_gauss_fermi_inf

! ***********************
! * INPUT variables     *
! ***********************
      real(dp), intent(in) :: E1, E2  ! energy parameters 
      real(dp), intent(in) :: kT      ! temperature in Ry
      integer,  intent(in) :: N

! ***********************
! * OUTPUT variables    *
! ***********************
      integer,  intent(out) :: N_mid, N_tail, GF_N_kT

! ***********************
! * LOCAL variables     *
! ***********************
      real(dp) :: deltaN, VkT

      VkT = dabs(E2-E1) / kT

      ! Determine the optimal Gauss-Fermi curve
      ! If the bias is larger than 10 * kT then we can select from [-5kT;inf]
      ! almost without overlap
      if ( VkT > 5._dp ) then
         ! we default to select -5kT (the fermi function is almost zero there)
         GF_N_kT = -5
      else
         ! We need to figure out were it makes sense to cut them
         ! Calculate the overlap of the fermi functions
         ! (5 - Vkt)/2 = overlap of Fermi functions
         ! add 1/2 to take the ceiling closest integer
         GF_N_kT = nint((5._dp-VkT)*.5_dp + .5_dp)
      end if

      ! Select an allowed range (this should not be a problem)
      GF_N_kT = min(GF_N_kT,G_NF_MAX_kT)
      GF_N_kT = max(GF_N_kT,G_NF_MIN_kT)

      ! Now we need to determine the number of points in each segment
      ! We do this by assuming that the weight in the tails
      ! equals that of the overlap region of the Fermi function.
      ! This seems like a reasonable assumption.

      ! Calculate the parts in the tail
      ! this is the delta E in the bias window
      deltaN = VkT / real(N,dp)

      if ( GF_N_kT > 0 ) then
         ! When extending the integral into the positive region
         ! it will be better to allow the full range to be divided
         deltaN = (VkT + GF_N_kT/2) / real(N,dp)
         N_tail = nint(.5_dp*GF_N_kT/deltaN)
      else
         N_tail = nint(real(abs(GF_N_kT),dp)/deltaN)
      end if

      ! correct for spill-out
      N_tail = min(N_tail,G_NF_MAX_N)
      N_tail = max(N_tail,G_NF_MIN_N)

      ! Calculate the number of points in the middle segment
      N_mid  = N - N_tail

    end subroutine init_GaussFermi

  end subroutine ts_read_contour_options


  subroutine ts_print_contour_options(cEq,cnEq,Eq_Eta,nEq_Eta,N_poles,IsVolt)

    use parallel, only : IONode
    use units, only : eV

    type(ts_io_c), intent(in) :: cEq(:), cnEq(:)
    real(dp), intent(in) :: Eq_Eta, nEq_Eta
    integer, intent(in) :: N_poles
    logical, intent(in) :: IsVolt

    character(len=200) :: chars
    character(len=200), parameter :: OPT_N = '(''ts_options: '',a)'
    character(len=200), parameter :: OPT_C = '(''ts_options: '',a,t53,''=    '',a)'
    character(len=200), parameter :: OPT_F = '(''ts_options: '',a,t53,''='',f10.4)'
    character(len=200), parameter :: OPT_INT = '(''ts_options: '',a,t53,''='',i5)'
    character(len=200), parameter :: OPT_F_U = '(''ts_options: '',a,t53,''='',f10.4,tr1,a)'
    character(len=200), parameter :: OPT_G_U = '(''ts_options: '',a,t53,''='',g11.4,tr1,a)'
    integer :: i

    if ( .not. IONode ) return
    
    write(*,opt_n) ' ----------------- Contour ----------------- '

    write(*,opt_n) '          >> Equilibrium contour << '
    write(*,opt_g_u) 'Equilibrium Greens function Eta',Eq_Eta/eV,'eV'
    do i = 1 , size(cEq)
       if ( i == 1 ) then
          chars = 'Circle'
       else if ( i == size(cEq) ) then
          chars = 'Line tail'
       else
          chars = 'Line contour'
       end if
       call write_e(trim(chars)//' contour E_min',cEq(i)%a)
       call write_e(trim(chars)//' contour E_max',cEq(i)%b)
       write(*,opt_int) trim(chars)//' contour points',cEq(i)%N
       write(*,opt_c) trim(chars)//' contour method', &
            trim(longtype2str(cEq(i)%type))
    end do
    write(*,opt_int) 'Number of poles', N_poles

    if ( IsVolt ) then
       write(*,opt_n) '        >> non-Equilibrium contour << '
       write(*,opt_g_u) 'non-Equilibrium Greens function Eta',nEq_Eta/eV,'eV'
       do i = 1 , size(cnEq)
          if ( i == 1 ) then
             chars = 'Lower tail'
          else if ( i == size(cnEq) ) then
             chars = 'Upper tail'
          else
             chars = 'Line'
          end if
          call write_e(trim(chars)//' contour E_min',cnEq(i)%a)
          call write_e(trim(chars)//' contour E_max',cnEq(i)%b)
          write(*,opt_int) trim(chars)//' contour points',cnEq(i)%N
          write(*,opt_c) trim(chars)//' contour method', &
               trim(longtype2str(cnEq(i)%type))
       end do

    end if

  contains

    subroutine write_e(str,val)
      character(len=*), intent(in) :: str
      real(dp), intent(in) :: val
      ! This is our definition of infinity....
      if ( abs(val) > 10000._dp ) then
         if ( val < 0._dp ) then
            write(*,opt_c) trim(str),'-Infinity'
         else
            write(*,opt_c) trim(str),' Infinity'
         end if
      else
         write(*,opt_f_u) trim(str),val / eV,'eV'
      end if
    end subroutine write_e
    
  end subroutine ts_print_contour_options

  subroutine ts_print_contour_warnings(cEq,cnEq,kT,Eq_Eta, nEq_Eta, N_poles, IsVolt)

    use parallel, only : IONode, Nodes, operator(.parcount.)
    use units, only : Pi, Kelvin

    type(ts_io_c), intent(in) :: cEq(:), cnEq(:)
    real(dp), intent(in) :: kT, Eq_Eta, nEq_Eta
    integer, intent(in) :: N_poles
    logical, intent(in) :: IsVolt
    integer :: i
    real(dp) :: diff

    if ( .not. IONode ) return

    call calc_Eta_Pole_Diff_Kelvin(Eq_Eta,diff)

    ! say we warn if we are 3 kT away from the pole
    if ( abs(diff) <= 10._dp ) then
       write(*,'(a)') 'NOTICE : Equilibrium Eta value is &
            &less than 10 Kelvin away from a pole.'
       write(0,'(a)') 'NOTICE : Equilibrium Eta value is &
            &less than 10 Kelvin away from a pole.'
    end if
    if ( abs(diff) <= 5._dp ) then
       write(*,'(a)') 'WARNING: Equilibrium Eta value is &
            &less than  5 Kelvin away from a pole.'
       write(0,'(a)') 'WARNING: Equilibrium Eta value is &
            &less than  5 Kelvin away from a pole.'
    end if

    if ( IsVolt ) then

       call calc_Eta_Pole_Diff_Kelvin(nEq_Eta,diff)

       ! say we warn if we are 3 kT away from the pole
       if ( abs(diff) <= 10._dp ) then
          write(*,'(a)') 'NOTICE : non-Equilibrium Eta value is &
               &less than 10 Kelvin away from a pole.'
          write(0,'(a)') 'NOTICE : non-Equilibrium Eta value is &
               &less than 10 Kelvin away from a pole.'
       end if
       if ( abs(diff) <= 5._dp ) then
          write(*,'(a)') 'WARNING: non-Equilibrium Eta value is &
               &less than  5 Kelvin away from a pole.'
          write(0,'(a)') 'WARNING: non-Equilibrium Eta value is &
               &less than  5 Kelvin away from a pole.'
       end if

       i = 2 * ( sum(cEq(:)%N) + N_poles )
       if ( mod(i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Equilibrium energy contour points are not"
          write(*,*) "        divisable by the number of nodes."
          write(*,*) "        Better scalability is achived by changing:"
          write(*,*) "          - TS.Contour.Eq.Circle.N"
          write(*,*) "          - TS.Contour.Eq.Line.N"
          write(*,*) "          - TS.Contour.Eq.Pole.N"
          write(*,*) "          - %block TS.Contour.Eq"

          ! Calculate optimal number of energy points
          write(*,'(t10,a,i4)') "Used equilibrium # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
               "Optimal equilibrium # of energy points: ",i, &
               "+- i*",Nodes
       end if
       
       i = sum(cnEq(:)%N)
       if ( mod(i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Non-equilibrium energy contour points are not"
          write(*,*) "        divisable by the number of nodes."
          write(*,*) "        Better scalability is achieved by changing:"
          write(*,*) "          - TS.Contour.nEq.N"
          write(*,*) "          - %block TS.Contour.nEq"
          
          ! Calculate optimal number of energy points
          write(*,'(t10,a,i4)') "Used non-equilibrium # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
               "Optimal non-equilibrium # of energy points: ",i, &
               "+- i*",Nodes
       end if
       
       i = 2 * sum(cEq(:)%N) + sum(cnEq(:)%N)
       if ( mod(i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Total energy contour points are not"
          write(*,*) "        divisable by the number of nodes."
          
          ! Calculate optimal number of energy points
          write(*,'(t10,a,i4)') "Used # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4,tr1,a4,i3,/)') &
               "Optimal # of energy points: ",i,"+- i*",Nodes
       end if
    else
       i = sum(cEq(:)%N) + N_poles
       
       !   - The equilibrium parts are the same computational cost
       !   * Solution make the equi contours divisible by Nodes
       if ( mod(i,Nodes) /= 0 ) then
          write(*,*) "NOTICE: Equilibrium energy contour points are not"
          write(*,*) "        divisable by the number of nodes."
          write(*,*) "        Better scalability is achived by changing:"
          write(*,*) "          - TS.Contour.Eq.Circle.N"
          write(*,*) "          - TS.Contour.Eq.Line.N"
          write(*,*) "          - TS.Contour.Eq.Pole.N"
          write(*,*) "          - %block TS.Contour.Eq"

          ! Calculate optimal number of energy points
          write(*,'(t10,a,i4)') "Used # of energy points   : ",i
          i = Nodes .PARCOUNT. i
          write(*,'(t10,a,i4)') "Optimal # of energy points: ",i
       end if

    end if

  contains

    subroutine calc_Eta_Pole_Diff_Kelvin(Eta,diff)
      real(dp), intent(in) :: Eta
      real(dp), intent(out) :: diff
      integer :: iPole

      ! We check whether the equilibrium shift into the imaginary plane
      ! is too close to a Fermi-pole
      ! We define this to be the case if 
      iPole = nint(Eta / (Pi * kT) )
      if ( iPole < 0 ) call die('Error in the eta value')
      if ( mod(iPole,2) == 0 ) then
         if ( iPole == 0 ) then
            ! we have the closest pole at kT * Pi
            diff = kT * Pi - Eta
         else
            ! we have the closest pole at either
            ! iPole - 1 or iPole + 1
            diff = Eta / ( kT * Pi ) 
            
            if ( abs(diff - real(iPole-1,dp)) < &
                 abs(diff - real(iPole+1,dp)) ) then
               diff = Eta - kT * Pi * real(iPole-1,dp)
            else
               diff = Eta - kT * Pi * real(iPole+1,dp)
            end if
         end if
         diff = abs(diff) / Kelvin
      else
         ! The eta value lies close to a pole
         ! check how close...
         diff = abs( Eta - kT * Pi * real(iPole,dp) ) / Kelvin
      end if

    end subroutine calc_Eta_Pole_Diff_Kelvin

  end subroutine ts_print_contour_warnings

  subroutine ts_read_contour_block(bName,cN,c, kT, V_half) 

    use parallel, only : IONode
    use units, only : eV
    use fdf

    character(len=*), intent(in) :: bName
    integer, intent(inout) :: cN
    type(ts_io_c), allocatable, intent(inout) :: c(:)
    real(dp), intent(in) :: kT
    real(dp), intent(in), optional :: V_half

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    character(len=200) :: g
    integer :: pNames, my_loop, cI, i, iV
    logical :: read_method, read_points, read_int, error

    ! if the block does not exist, return
    if ( .not. fdf_block(bName,bfdf) ) return

    ! counter for the number of lines with "points"
    cN = 0

    ! first read in the amount of points we need...
    do while ( fdf_bline(bfdf,pline) )
       
       ! figure out if we have points in the line
       if ( fdf_bnnames(pline) > 0 ) then
          
          ! Read the first name in the line
          g = fdf_bnames(pline,1)
          
          if ( leqi(g,'points') .or. &
               leqi(g,'pts') .or. &
               leqi(g,'n') ) then
             cN = cN + 1
          end if

       end if

    end do

    if ( cN == 0 ) &
         call die('No contour segments found in block: ' &
         //trim(bName)//'. Please correct input')

    ! We now have the number of points in the block
    if ( allocated(c) ) deallocate(c)
    allocate(c(cN))

    ! move back to be read again
    call fdf_brewind(bfdf)

    ! Read the information in block

    ! cI is the current contour segment.
    cI = 0
    rCB: do while ( fdf_bline(bfdf,pline) )

       ! if no names exist we can loop to the next line
       ! thus commenting out a contour works very good!
       if ( fdf_bnnames(pline) == 0 ) cycle rCB

       ! increment the contour
       cI = cI + 1

       ! ensure that the contour type is not defined
       ! important for the Gauss-Fermi check
       c(cI)%type = 0
       ! default the optional argument to not being set
       c(cI)%opt     = huge(0._dp)
       c(cI)%G_F_off = 0._dp
       
       read_method   = .false.
       read_points   = .false.
       read_int      = .false.

       ! we have a three line segment now...
       ! we allow the user to do it as (s)he pleases...
       read_CS: do my_loop = 1 , 3

       ! we read the next line in loop 2 and 3
       if ( my_loop > 1 ) then
          if ( cI < cN .and. .not. fdf_bline(bfdf,pline) ) then
             call die('Block: '//trim(bName)//' is not build correctly.')
          end if
       end if

       pNames = fdf_bnnames(pline)
       ! figure out if we have any names in the line
       ! all 3 lines *must* have at least one name in the line.
       if ( pNames == 0 ) &
            call die('No names are present in the line &
            &please correct the block: '//trim(bName))

       i = 0
       do while ( i <= pNames ) 
          
          i = i + 1

          ! Read the name in the line
          g = fdf_bnames(pline,i)

          if ( leqi(g,'points') .or. &
               leqi(g,'pts') .or. &
               leqi(g,'n') ) then

             if ( fdf_bnintegers(pline) < 1 ) then
                call die('Block '//trim(bName)//' is erroneous &
                     &could not figure out the number of points...')
             end if

             if ( read_points ) call die('Block '//trim(bName)//' is erroneous. &
                  &Two point specifications are found for one contour segment.')

             ! the number of points on that contour is always the first
             ! integer
             c(cI)%N = fdf_bintegers(pline,1)

             ! ensure that we have read the number of points segment for 
             ! the current contour segment...
             read_points = .true.

             cycle read_CS

          else if ( leqi(g,'method') .or. &
               leqi(g,'m') ) then

             if ( read_method ) call die('Block '//trim(bName)//' is erroneous. &
                  &Two methods are found for one contour segment.')
             
             if ( i == pNames ) call die('The method is not specified in block &
                  &'//trim(bName)//'. Please correct your FDF block.')

             c(cI)%type = str2type(fdf_bnames(pline,i+1))

             ! if any values exist we save it... (it is an optional argument)
             if ( fdf_bnvalues(pline) > 0 ) then
                c(cI)%opt = fdf_bvalues(pline,1)
             end if
             
             if ( c(cI)%type == CC_TYPE_G_NF_0kT ) then
                ! check that the Gauss-fermi contour is the first or
                ! last segment
                if ( 1 < cI .and. cI < cN ) then
                   call die('Gauss-Fermi contour segment can only be placed &
                        &in the tails. Please correct block: '//trim(bName))
                end if

             end if

             ! We need to correct for the Gauss-Fermi contour...
             call correct_G_Fermi(c(cI),kT)

             ! If the previous one is "next" we can set it
             if ( 1 < cI ) then
                if ( leqi(c(cI-1)%c_b,'next') ) then
                   c(cI-1)%b = c(cI)%a
                end if
             end if

             ! ensure that we have read the number of points segment for 
             ! the current contour segment...
             read_method = .true.

             cycle read_CS

          else if ( leqi(g,'to') .or. &
               leqi(g,'--') .or. &
               leqi(g,'---') ) then

             if ( read_int ) call die('Block '//trim(bName)//' is erroneous. &
                  &Two integration bounds are found for one contour segment.')

             ! This is the "tricky" setting
             ! In case the first "number" is -inf we need to specify the correct 
             ! value number on the line...
             iV = 1
             ! The user can specify previous and next

             ! first we check whether previous is set...
             if ( i == 1 ) then
                call die('The integrate line must have at least a unit or &
                     &another keyword: inf|previous as the first input. &
                     &Error in block: '//trim(bName))
             else if ( i == pNames ) then
                call die('The integrate line must have at least a unit or &
                     &another keyword: inf|next as the last input. &
                     &Error in block: '//trim(bName))
             end if

             ! this should be one of: next|inf|<unit>
             g = fdf_bnames(pline,i-1)
             if ( leqi(g,'previous') .or. leqi(g,'prev') ) then
                if ( cI < 2 ) call die('You cannot specify previous &
                     &for the first contour. Please correct block: '//trim(bName))

                ! Set the previous value correctly...
                c(cI)%c_a = 'prev'
                c(cI)%a = c(cI-1)%b
                
             else

                ! lower integration bounds has not been set
                ! figure out what it is...

                ! We need to parse the format...
                call pline_E_parse(pline,1,i-1,iV, &
                     c(cI), .true., &
                     V_half=V_half)

             end if

             ! this should be one of: next|inf|<unit>
             g = fdf_bnames(pline,i+1)
             if ( leqi(g,'next') ) then
                if ( cI == cN ) call die('You cannot specify next &
                     &for the last contour segment. &
                     &Please correct block: '//trim(bName))

                c(cI)%c_b = 'next'

             else

                ! upper integration bounds has not been set
                ! figure out what it is...
                call pline_E_parse(pline,i+1,pNames,iV, &
                     c(cI), .false., &
                     V_half=V_half)
             end if

             ! we need to correct for the Gauss-Fermi method... (in case method
             ! was set before)
             call correct_G_Fermi(c(cI), kT)

             ! Set the previous one is "next" we can set it
             if ( 1 < cI ) then
                if ( leqi(c(cI-1)%c_b,'next') ) then
                   c(cI-1)%b = c(cI)%a
                end if
             end if

             read_int = .true.

             cycle read_CS

          end if

       end do

       end do read_CS

       if ( .not. ( read_int .and. read_method &
            .and. read_points ) ) then
          call die('Have not succesfully read a full contour segment &
               &in block '//trim(bName)//' please correct the block.')
       end if

       
    end do rCB

    if ( cI /= cN ) &
         call die('Unsuccesfully read the amount of segments that was &
         &initially found. Please correct block: '//trim(bName))

  contains 

    ! This subroutine corrects and checks the Gauss-Fermi input
    ! we restrict the user to be accurate for the energy designation
    ! of the kT separation from the V/2 points.
    ! If we don't find any integer valued separation from the common
    ! mean, then the job is killed.
    ! Lastly we set the c%opt to be the bottom-end integral
    subroutine correct_G_Fermi(c, kT)
      type(ts_io_c), intent(inout) :: c
      real(dp)  , intent(in)    :: kT

      real(dp), parameter :: kT_error = 1.e-6_dp
      real(dp) :: r_a,r_b
      logical :: error
      integer :: a_i, b_i

      ! if not a Gauss-Fermi contour, return...
      if ( c%type /= CC_TYPE_G_NF_0kT ) return

      a_i = huge(1)
      b_i = huge(1)

      ! First correct the integer kT resolution.
      ! we do not allow for other kT resolutions...

      ! an error tolerance of 1e-6 Ry should be "very fine"
      ! it is on the order of 1e-7 eV
      error = .false.

      r_a = c%a
      r_b = c%b

      ! We must ensure that the infinite points
      ! are not gobbled...
      ! we must not create Gauss-Fermi quadrature points for \int_x^10000kT
      ! but why should we ever need that? :)
      if ( abs(c%a) / kT < 10000._dp ) then
         a_i = nint((c%a - c%G_F_off)/kT)
         c%a = kT * a_i + c%G_F_off
         error = error .or. ( abs(c%a - r_a) > kT_error )
      end if
      if ( abs(c%b) / kT < 10000._dp ) then
         b_i = nint((c%b - c%G_F_off)/kT)
         c%b = kT * b_i + c%G_F_off
         error = error .or. ( abs(c%b - r_b) > kT_error )
      end if

      ! Correct for the off-set
      r_a = r_a - c%G_F_off
      r_b = r_b - c%G_F_off

      if ( error ) then
         write(*,'(a25,'': '',g20.13)') 'Offset',c%G_F_off / eV
         write(*,'(a25,'': '',a)'     ) 'Lower string formatting',trim(c%c_a)
         write(*,'(a25,'': '',g20.13)') 'Original lower bound',r_a / eV
         write(*,'(a25,'': '',g20.13)') 'Truncated lower bound',(c%a - c%G_F_off)/ eV
         write(*,'(a25,'': '',i0)'    ) 'Lower kT integer found',a_i
         write(*,'(a25,'': '',g20.13)') 'kT contribution',a_i * kT
         write(*,*)
         write(*,'(a25,'': '',a)'     ) 'Upper string formatting',trim(c%c_b)
         write(*,'(a25,'': '',g20.13)') 'Original upper bound',r_b/eV
         write(*,'(a25,'': '',g20.13)') 'Truncated upper bound',(c%b - c%G_F_off)/eV
         write(*,'(a25,'': '',i0)'    ) 'Upper kT integer found',b_i
         write(*,'(a25,'': '',g20.13)') 'kT contribution',b_i * kT

         write(*,'(2(/,a))') 'Energy specification for the Gauss-Fermi &
              &function does not fully correspond to an integer kT separation.', &
              'Please only use V/2 and <int> kT or inf for the contour of the &
              &Gauss-Fermi functions.'
         call die('')
      end if
      
      ! We need to ensure that the end points are valid...
      ! first determine which of the ends is the largest...
      if ( abs(c%a) < abs(c%b) ) then
         ! We have an integral in the positive end

         c%type = CC_TYPE_G_NF_0kT + a_i
         c%G_F_lower = a_i
         c%G_F_upper = huge(1)
         if ( b_i < huge(1) ) then
            c%G_F_upper = b_i
         end if

         ! We might as well ensure that a_i is correct here...
         select case ( c%G_F_lower ) 
         case ( G_NF_MIN_kT : G_NF_MAX_kT )
               ! do nothing, valid input
         case default
            call die('The requested Gauss-Fermi contour does &
                 &not exist. Please correct the block: '//trim(bName))
         end select
         
      else
         ! We have an integral in the negative end

         c%type = CC_TYPE_G_NF_0kT - b_i
         c%G_F_lower = -b_i
         c%G_F_upper = huge(1)
         if ( abs(a_i) < huge(1) ) then
            c%G_F_upper = abs(a_i)
         end if

         ! We might as well ensure that a_i is correct here...
         select case ( c%G_F_lower ) 
         case ( G_NF_MIN_kT : G_NF_MAX_kT )
               ! do nothing, valid input
         case default
            call die('The requested Gauss-Fermi contour does &
                 &not exist. Please correct the block: '//trim(bName))
         end select
         
      end if

    end subroutine correct_G_Fermi

    subroutine pline_E_parse(pline,iS,iE,iV,cc,is_a,V_half)
      ! This routine parses the pline and extracts the energy
      ! which is derived from the expression...
      type(parsed_line), pointer :: pline
      integer, intent(in) :: iS, iE
      integer, intent(inout) :: iV
      type(ts_io_c), intent(inout) :: cc
      logical, intent(in) :: is_a
      real(dp), intent(in), optional :: V_half

      real(dp) :: E, adjust
      character(len=c_ab_N) :: c
      character(len=50) :: g
      logical :: add
      integer :: i, addV, cL

      ! initialize parsing...
      add = .true.
      E = 0._dp
      addV = iV
      adjust = 0._dp
      c = ' '

      ! begin parsing...
      do i = iS, iE

         g = fdf_bnames(pline,i)

         ! Parse single elements
         if ( leqi(g,'+') ) then
            add = .true.
            cycle

         else if ( leqi(g,'-') ) then
            add = .false.
            cycle

         end if

         ! Parse compressed elements
         if ( leqi(g(1:1),'+') ) then
            add = .true.
            g = g(2:)

         else if ( leqi(g(1:1),'-') ) then
            add = .false.
            g = g(2:)

         end if

         ! we need to be sure to capture everything...
         if ( leqi(g,'inf') ) then
            if ( add ) then
               E = huge(0._dp)
               c = 'inf'
            else
               E = -huge(0._dp)
               c = '-inf'
            end if
         else if ( leqi(g,'kt') .or. leqi(g,'kbt') ) then
            if ( add ) then
               E = E + fdf_bvalues(pline,iV) * kT
               write(c,'(a,'' + '',g20.13,tr1,a)') &
                    trim(c),fdf_bvalues(pline,iV), 'kT'
            else
               E = E - fdf_bvalues(pline,iV) * kT
               write(c,'(a,'' - '',g20.13,tr1,a)') &
                    trim(c),fdf_bvalues(pline,iV), 'kT'
            end if
            iV = iV + 1
         else if ( leqi(g,'v/2') ) then
            if ( .not. present(V_half) ) then
               call die('Requesting the voltage is not supported in the &
                    &block: '//trim(bName))
            end if
            if ( abs(adjust) /= 0._dp ) then
               call die('You can only use one V_half reference per &
                    &energy.')
            end if

            if ( add ) then
               E = E + abs(V_half)
               adjust =   abs(V_half)
               write(c,'(a,'' + V/2'')') trim(c)
            else
               E = E - abs(V_half)
               adjust = - abs(V_half)
               write(c,'(a,'' - V/2'')') trim(c)
            end if

         else 
            ! We neet to check that 'g' is actually an energy unit...
            if ( .not. ( leqi(g,'ev') .or. leqi(g,'mev') .or. &
                 leqi(g,'ry') .or. leqi(g,'mry') ) ) then
               ! It is not a proper energy unit...
               ! hence the value most not exist...
               ! TODO add an fdf-check for the location of the next
               ! value. If it is preceding this name we should error
               ! out in that the unit is wrong!
               ! However, currently the printing of the block
               ! shows the difference...
               cycle
            end if

            if ( add ) then
               E = E + fdf_bvalues(pline,iV) * fdf_convfac(g,'Ry')
               write(c,'(a,'' + '',g20.13,tr1,a)') &
                    trim(c),fdf_bvalues(pline,iV), trim(g)
            else
               E = E - fdf_bvalues(pline,iV) * fdf_convfac(g,'Ry')
               write(c,'(a,'' - '',g20.13,tr1,a)') &
                    trim(c),fdf_bvalues(pline,iV), trim(g)
            end if
            iV = iV + 1
         end if

         ! in case the string gets bloated with a lot of options
         call correct_E_string(c)

         ! We default to add
         ! in case the sign is part of the value we need to add
         add = .true.

      end do

      if ( abs(adjust) > 0._dp .and. iV - addV > 1 ) then
         call die('You can only use the voltage specification and &
              &one other expression. Please simplify your input.')
      end if

      if ( cc%G_F_off == 0._dp ) cc%G_F_off = adjust
      if ( is_a ) then
         cc%a = E
         cc%c_a = c
      else
         cc%b = E
         cc%c_b = c
      end if

    end subroutine pline_E_parse

  end subroutine ts_read_contour_block

  subroutine correct_E_string(c)
    character(len=*), intent(inout) :: c
    integer :: i, cL

    ! We need to correct the 0'es in the character line...
    cL = len_trim(c)
    do i = cL - 1 , 1 , -1
       ! trim '0 ' and double spaces '  '
       if ( c(i:i+1) == '0 ' .or. c(i:i+1) == '  ' ) then
          c(i:cL-1) = c(i+1:cL)
          cL = cL - 1
       end if

       ! trim double signs... '+ -' or '- +'
       if ( i < cL - 2 ) then
          if ( c(i:i+2) == '+ -' ) then
             c(i:i) = '-'
             c(i+2:cL-1) = c(i+3:cL)
             cL = cL - 1
          else if ( c(i:i+2) == '- +' ) then
             c(i+2:cL-1) = c(i+3:cL)
             cL = cL - 1
          end if
       end if
    end do

    ! trim leading ' ' and '+'
    do while ( c(1:1) == ' ' .or. c(1:1) == '+' )
       c(1:cL-1) = c(2:cL)
       cL = cL - 1
    end do

    c(cL+1:) = ' '

  end subroutine correct_E_string

  subroutine ts_print_contour_block(bName,c,msg_first,msg_last)

    use parallel, only : IONode

    character(len=*), intent(in) :: bName
    type(ts_io_c), intent(in) :: c(:)
    character(len=*), intent(in), optional :: msg_first, msg_last
    integer :: i, cN

    if ( .not. IONode ) return

    ! Start by writing out the block beginning
    write(*,'(a,a)') '%block ',trim(bName)

    ! Now we write out the information
    cN = size(c)
    do i = 1 , cN
       if ( 1 < i ) write(*,'(a)') ''
       if ( i == 1 .and. present(msg_first) ) then
          write(*,'(t3,a)') trim(msg_first)
       else if ( i == cN .and. present(msg_last) ) then
          write(*,'(t3,a)') msg_last
       end if

       call write_contour_block_E_part(c(i))
       
    end do
    
    write(*,'(a,a)') '%endblock ',trim(bName)
    
  contains 
    
    subroutine write_contour_block_E_part(c)
      type(ts_io_c), intent(in) :: c

      ! move to third column...
      write(*,'(t3)',advance='no')
      
      write(*,'(3a)',advance='no') trim(c%c_a),' to ',trim(c%c_b)
      
      write(*,'(a)') ' with'
      
      ! Print the number of points...
      write(*,'(t7,a,tr1,i0,tr1,a)') 'points', c%N, 'using'
      
      ! Print the method
      write(*,'(t9,a,tr1,a)') 'method', trim(type2input(c%type))
      
    end subroutine write_contour_block_E_part
    
  end subroutine ts_print_contour_block


  ! Routine for "pretty" printing the contour points out in the out file
  subroutine ts_print_contour(contour)
    use parallel, only : IONode
    use units,    only : eV

    type(ts_ccontour), pointer :: contour(:)
! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c
    character(len=CC_TYPE_LEN) :: ctype
    character(len=3) :: f_type
    integer :: i, part, type

    if ( .not. IONode ) return
    
    ! Initialize variables
    part = -1
    type = -1
    nullify(c)

    write(*,'(a)') "transiesta: contour integration path:"
    write(f_type,'(a,i2)') 'a',CC_TYPE_LEN
    write(*,'(1x,'//trim(f_type)//',''   '',2(tr1,a12),2(tr1,a14))') &
         "Type  ","Re(c)[eV]","Im(c)[eV]","Re(weight)","Im(weight)"
    do i = 1 , size(contour)
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

  end subroutine ts_print_contour


! Write out the contour to a contour file
  subroutine ts_io_contour(contour,slabel,suffix)
    use parallel, only : IONode
    use units, only : eV
    type(ts_ccontour), intent(in) :: contour(:)
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

! *********************
! * LOCAL variables   *
! *********************
    character(len=len_trim(slabel)+9) :: fname
    integer :: i, unit, part

    if ( .not. IONode ) return

    if ( present(suffix) ) then
       fname = trim(trim(slabel)//trim(suffix))
    else
       fname = trim(trim(slabel)//'.TSCC')
    end if
    call io_assign( unit )
    open( unit, file=fname, status='unknown' )
    write(unit,'(a)') '# Complex contour path'
    write(unit,'(a,a12,3(tr1,a13))') '#','Re(c)[eV]','Im(c)[eV]','Re(w)','Im(w)'
    part = contour(1)%part
    do i = 1 , size(contour)
       if ( part /= contour(i)%part ) then
          write(unit,*)
          part = contour(i)%part
       end if
       write(unit,'(4(e13.6,tr1))') contour(i)%c/eV,contour(i)%w
    end do
    call io_close( unit )

  end subroutine ts_io_contour

end module m_ts_io_contour
