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
!   6) transport
!   7) phonon

! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_cctype
  use m_ts_io_contour
  use precision, only : dp

  implicit none

  ! This module will also contain all the contour variables
  integer, save :: NEn_eq
  integer, save :: NEn  ! Number of energy points in the contour path
  integer, save :: PNEn ! Number of energy points in the contour path (divisible by Nodes)

  ! Contour path
  type(ts_ccontour), dimension(:), pointer, save :: contour
  type(ts_io_c), allocatable, public :: cEq(:), cnEq(:)

  ! The contour specific variables
  real(dp), save, public :: kT, EfL, EfR
  real(dp), save, public :: nEq_Eta, Eq_Eta

  ! We need to retain the information about the contour here.
  ! It provides an easier overview as there are quite a few constants governing the
  ! methods.
  integer, save, public  :: N_poles, N_transport

  ! The contours for the equilibrium density are attributed a fruitful discussion with
  ! Hans Skriver. Previously the routine names reflected his contribution.
  ! However, for programming clarity we have employed a different naming scheme.
  ! In order to retain the contributions it is encouraged to keep this sentiment for his
  ! contributions.

  ! Furthermore the non-equilibrium density integration are attributed discussions with
  ! Antti-Pekka Jauho. 

  ! For further attributions see the original paper by Brandbyge, et. al, 2002: DOI: 10.1103/PhysRevB.65.165401

  public :: NEn, PNEn, contour
  public :: setup_contour
  public :: sort_contour
  public :: contour_Eq
  public :: contour_EqL, contour_EqR, contour_nEq
  public :: contour_Transport
  public :: read_contour_options

  private

contains

  subroutine read_contour_options(kT_in, IsVolt, VoltL, VoltR)

    use units, only : eV
    use parallel, only : IONode, Nodes, operator(.parcount.)
    use fdf

    logical, intent(in) :: IsVolt
    real(dp), intent(in) :: kT_in, VoltL, VoltR ! in Ry
    logical :: correct
    integer :: i,j

    ! Transfer options for the fermi-levels
    EfL = VoltL
    EfR = VoltR

    ! Copy over temperature
    kT = kT_in

    ! Read in the generic things about the contours...
    call fdf_deprecated('TS.ComplexContour.NPoles','TS.Contour.Eq.Pole.N')
    N_poles = fdf_get('TS.ComplexContour.NPoles',6)
    N_poles = fdf_get('TS.Contour.Eq.Pole.N',N_poles)
    
    ! broadening
    Eq_Eta = fdf_get('TS.Contour.Eq.Eta',0._dp,'Ry')
    if ( Eq_Eta < 0._dp) call die('ERROR: Eq_Eta < 0, we do not allow &
         &for using the advanced Greens function, please correct.')

    call fdf_deprecated('TS.biasContour.Eta','TS.Contour.nEq.Eta')
    nEq_Eta = fdf_get('TS.biasContour.Eta',0.000001_dp*eV,'Ry')
    nEq_Eta = fdf_get('TS.Contour.nEq.Eta',nEq_Eta,'Ry')
    if ( nEq_Eta <= 0._dp ) call die('ERROR: nEq_Eta <= 0, we do not allow &
         &for using the advanced Greens function, please correct.')

    ! Read in everything about the contour.
    call ts_read_contour_options(cEq,cnEq,kT,IsVolt, V_half=VoltL)

    N_transport = 0

    ! Correct for imperfect computational balance
    correct = fdf_get('TS.Contour.Eq.NoEmptyCycles',.true.)
    correct = correct .and. ( mod(sum(cEq(:)%N),Nodes) == 0 )
    if ( correct ) then
       i = sum(cEq(:)%N) + N_poles
       if ( IsVolt ) i = i * 2
       j = 1
       do while ( mod(i,Nodes) /= 0 )
          cEq(j)%N = cEq(j)%N + 1
          i = i + 1
          j = j + 1
          if ( j >= size(cEq) - 1 ) j = 1
       end do
    end if

    correct = fdf_get('TS.Contour.nEq.NoEmptyCycles',.true.)
    if ( correct ) then
       i = sum(cnEq(:)%N)
       if ( size(cnEq) == 1 ) then
          ! it must be a Sommerfeld contour...
          cnEq(1)%N = cnEq(1)%N + Nodes - mod(i,Nodes)
       else
          j = 2
          do while ( mod(i,Nodes) /= 0 )
             cnEq(j)%N = cnEq(j)%N + 1
             i = i + 1
             j = j + 1
             if ( j >= size(cnEq) ) j = 2
          end do
       end if
    end if

    ! Update the number of contour points
    NEn_eq = sum(cEq(:)%N) + N_poles
    NEn = NEn_eq 
    if ( IsVolt ) then
       NEn = NEn_eq * 2
       NEn = NEn + sum(cnEq(:)%N)
    end if
    
    PNEn = Nodes .parcount. NEn

  end subroutine read_contour_options

  ! Routine for creating the contour
  subroutine setup_contour(IsVolt)
    
    use precision, only : dp
    use parallel,  only : IONode, Nodes, operator(.PARCOUNT.)
    use m_ts_aux, only : nf2

    use m_gauss_quad
    use m_gauss_fermi_inf
    use m_gauss_fermi_30
    use m_gauss_fermi_28
    use m_gauss_fermi_26
    use m_gauss_fermi_24
    use m_gauss_fermi_22
    use m_gauss_fermi_20
    use m_gauss_fermi_19
    use m_gauss_fermi_18
    use m_gauss_fermi_17

! **********************
! * INPUT variables    *
! **********************
    logical,  intent(in) :: IsVolt ! Do we have a volt

! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c(:) => null()
    real(dp), allocatable :: x(:), w(:)
    real(dp) :: sE1, sE2, tmp
    logical :: switched
    integer :: Net, N
    integer :: i, j, k, l

    ! ensure that the ar
    nullify(contour)
    allocate(contour(NEn))

    ! Setup all the different methods
    if ( .not. IsVolt ) then

       c => contour_Eq()
       ! We will only create the equilibrium contour
       call setup_contour_Eq(CC_PART_EQUI_CIRCLE, cEq(:), &
            N_poles, &
            0._dp, kT, Eq_Eta, c)

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
       call setup_contour_Eq(CC_PART_L_EQUI_CIRCLE, cEq(:), &
            N_poles, &
            EfL, kT, Eq_Eta, c)

       c => contour_EqR()
       ! We will only create the equilibrium contour
       call setup_contour_Eq(CC_PART_R_EQUI_CIRCLE, cEq(:), &
            N_poles, &
            EfR, kT, Eq_Eta, c)

       ! The last contour is the non-equilibrium
       
       ! The voltage contour
       c => contour_nEq()
       j = 0

       N = 0
       ! find the maximum number of points in a single contour segment.
       do i = 1 , size(cnEq)
          N = max(N,cnEq(i)%N)
       end do
       allocate(x(N),w(N))

       ! loop over all contour segments
       do k = 1 , size(cnEq(:))

          N = cnEq(k)%N

          if ( cnEq(k)%type == CC_TYPE_SOMMERFELD ) then
             ! The only special type of contour
             ! this one will not allow other contour segments
             
             call Sommerfeld(EfL, EfR, kT, N, x, w)
             
             if ( size(cnEq) > 1 ) then
                call die('You cannot have several contour &
                     &segments and use the Sommerfeld contour.')
             end if

          else
          
             select case ( cnEq(k)%type ) 
             case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )

                ! Determine the method
                select case ( cnEq(k)%G_F_upper )
                case ( huge(1) )
                   call GaussFermi_inf(cnEq(k)%G_F_lower,N,x,w)
                case ( 30 )
                   call GaussFermi_30(cnEq(k)%G_F_lower,N,x,w)
                case ( 28 )
                   call GaussFermi_28(cnEq(k)%G_F_lower,N,x,w)
                case ( 26 ) 
                   call GaussFermi_26(cnEq(k)%G_F_lower,N,x,w)
                case ( 24 ) 
                   call GaussFermi_24(cnEq(k)%G_F_lower,N,x,w)
                case ( 22 )
                   call GaussFermi_22(cnEq(k)%G_F_lower,N,x,w)
                case ( 20 )
                   call GaussFermi_20(cnEq(k)%G_F_lower,N,x,w)
                case ( 19 )
                   call GaussFermi_19(cnEq(k)%G_F_lower,N,x,w)
                case ( 18 )
                   call GaussFermi_18(cnEq(k)%G_F_lower,N,x,w)
                case ( 17 )
                   call GaussFermi_17(cnEq(k)%G_F_lower,N,x,w)
                case default
                   call die('Unknown tail integral ending')
                end select

                ! shift to the correct placement of the contour.
                if ( abs(cnEq(k)%a) > abs(cnEq(k)%b) ) then
                   do i = 1 , N / 2
                      l = N + 1 - i
                      tmp = x(i)
                      x(i) = x(l)
                      x(l) = tmp
                      tmp = w(i)
                      w(i) = w(l)
                      w(l) = tmp
                   end do
                end if

                do i = 1 , N
                   ! reverse sign in case of negative Gauss-Fermi
                   if ( abs(cnEq(k)%a) > abs(cnEq(k)%b) ) x(i) = -x(i)
                   x(i) = x(i) * kT + cnEq(k)%G_F_off
                   w(i) = w(i) * kT
                end do

             case default
                
                call line_integral(cnEq(k)%type,cnEq(k)%a, cnEq(k)%b, kT, N, x, w)
                
                do i = 1 , N
                   w(i) = w(i) * nf2(x(i),sE1,sE2,kT)
                end do
                
             end select
             
          end if

          ! We need to reverse the arguments here (positive bias is "standard")
          call correct_weight_sign(N,w,EfR,EfL)

          do i = 1 , N
          
             j = j + 1
             c(j)%c = dcmplx(x(i),nEq_Eta)
             c(j)%w = dcmplx(w(i),0._dp)
             
             c(j)%part = CC_PART_NON_EQUI
             c(j)%type = cnEq(k)%type
          end do

       end do

       deallocate(x,w)

    end if

    ! Currently the contour segment does not exist in transiesta
    ! we could add a transmission function to track the transmission.

  end subroutine setup_contour


  ! Routine for creating the contour
  subroutine setup_contour_Eq(PART, cEq, N_poles, Ef, kT, Eta, contour)

    use units, only : Pi, eV
    use precision, only : dp
    use parallel,  only : IONode, Nodes, operator(.PARCOUNT.)
    use m_ts_aux, only : nf

    use m_gauss_quad
    use m_gauss_fermi_inf
    use m_gauss_fermi_30
    use m_gauss_fermi_28
    use m_gauss_fermi_26
    use m_gauss_fermi_24
    use m_gauss_fermi_22
    use m_gauss_fermi_20
    use m_gauss_fermi_19
    use m_gauss_fermi_18
    use m_gauss_fermi_17

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: PART
    type(ts_io_c), intent(in) :: cEq(:)
    integer, intent(in) :: N_poles
    real(dp), intent(in) :: Ef, kT, Eta

    type(ts_ccontour), pointer, intent(out) :: contour(:)

! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c(:) => null()
    real(dp), allocatable :: x(:), w(:)
    real(dp) :: Delta, D, gamma, alpha, R, beta
    complex(dp) :: z0, ztmp
    integer :: N_kT, G_nf_END
    integer :: i, j, k, N, curN

    if ( size(contour) <= 0 ) then  ! if there are not any points we definetely have an error
       call die("ERROR: setup_contour_Eq: no energy points specified")
    else if ( size(contour) /= sum(cEq(:)%N) + N_poles ) then
       call die("ERROR: setup_contour_Eq: no energy points erroneous.")
    end if
    
    i = 1
    curN = N_poles
    if ( curN > size(contour) ) call die('Error in contour setup')
    ! Initialize the poles
    c => contour(1:curN)
    call nf_poles(PART+2,Ef,kT,Eta,N_poles,c)

    ! Add the line integral (i.e. from above the last pole)
    Delta = aimag(c(curN)%c) + Pi*kT
    
    ! Add all line integrals...
    do k = size(cEq) , 2 , -1
       N = cEq(k)%N
       j = curN + 1
       curN = curN + N
       if ( curN > size(contour) ) call die('Error in contour setup')
       if ( N <= 0 ) cycle
       c => contour(j:curN)

       ! Allocate space
       allocate(x(N),w(N))

       ! Determine the method
       select case ( cEq(k)%type ) 
       case ( CC_TYPE_G_NF_MIN:CC_TYPE_G_NF_MAX )

          select case ( cEq(k)%G_F_upper )
          case ( huge(1) )
             call GaussFermi_inf(cEq(k)%G_F_lower,N,x,w)
          case ( 30 )
             call GaussFermi_30(cEq(k)%G_F_lower,N,x,w)
          case ( 28 )
             call GaussFermi_28(cEq(k)%G_F_lower,N,x,w)
          case ( 26 ) 
             call GaussFermi_26(cEq(k)%G_F_lower,N,x,w)
          case ( 24 ) 
             call GaussFermi_24(cEq(k)%G_F_lower,N,x,w)
          case ( 22 )
             call GaussFermi_22(cEq(k)%G_F_lower,N,x,w)
          case ( 20 )
             call GaussFermi_20(cEq(k)%G_F_lower,N,x,w)
          case ( 19 )
             call GaussFermi_19(cEq(k)%G_F_lower,N,x,w)
          case ( 18 )
             call GaussFermi_18(cEq(k)%G_F_lower,N,x,w)
          case ( 17 )
             call GaussFermi_17(cEq(k)%G_F_lower,N,x,w)
          case default
             call die('Unknown tail integral ending')
          end select

       case default 
          call die('Nothing could be deciphered from the Eq. tail contour method')
       end select
    
       ! We move them over to this
       do i = 1 , N
          j = N + 1 - i
          c(i)%c    = dcmplx(x(j)*kT + Ef,Delta)
          c(i)%w    = - w(j) * kT
          c(i)%part = PART + 1
          c(i)%type = cEq(k)%type
       end do
       
       deallocate(x,w)
    end do

    ! We create the circle contour
    ! First we lift the equilibrium contour from the 
    ! real axis by Eta
    Delta = Delta - Eta
    ! The ending point for the contour (on the imaginary axis)
    gamma = cEq(1)%b
    D     = -cEq(1)%a + gamma 
    ! The angle between CCEmin and the point
    alpha = dATAN(Delta/D)
    ! The radius of the circle in the complex plane
    R     = dsqrt(Delta*Delta + D*D) / (2._dp*cos(alpha))
    ! The center of the circle on the real-axis shifted by Eta
    z0    = dcmplx(cEq(1)%a + Ef + R, Eta)
    ! The angle from where we start
    beta  = dasin(Delta / R)

    ! We now create the circle contour
    i = curN + 1
    curN = curN + cEq(1)%N
    if ( curN > size(contour) ) call die('Error in contour setup')
    c => contour(i:curN)
    
    ! Allocate space
    allocate(x(cEq(1)%N),w(cEq(1)%N))

    select case ( cEq(1)%type )
    case ( CC_TYPE_G_LEGENDRE )
       call Gauss_Legendre_Rec(cEq(1)%N, 0, beta, Pi, x, w)
    case ( CC_TYPE_TANH_SINH )
       ! TODO add the optional argument
       !D = 1.2e-2_dp * R / real(cEq(1)%N,dp)
       D = 1.8e-2_dp * (Pi-beta) / real(cEq(1)%N,dp)
       call TanhSinh_Exact(cEq(1)%N, x, w, beta, Pi, p=D)
    case default
       call die('Unrecognized contour for the equilibrium circle')
    end select

    do i = 1 , cEq(1)%N
       ztmp   = R * cdexp(dcmplx(0._dp,x(i)))
       c(i)%c = z0 + ztmp
       ! Factor i, comes from Ed\theta=dE=iR e^{i\theta}
       c(i)%w = w(i) * nf((c(i)%c-Ef)/kT) * dcmplx(0._dp,1._dp)*ztmp
       c(i)%part = PART
       c(i)%type = cEq(1)%type
    end do

    deallocate(x,w)

    ! Note we put a minus here because the integral we want is the
    ! negative of the pole-sum and C+L integral
    contour(:)%w = - contour(:)%w
    
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
    use m_integrate

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
    case ( CC_TYPE_SIMP_MIX )
      
       call Simpson_38_3_rule(NEn,x,w,E1,E2)
       
    case ( CC_TYPE_BOOLE_MIX )
      
       call Booles_Simpson_38_3_rule(NEn,x,w,E1,E2)
       
    case ( CC_TYPE_MID )

       call Mid_Rule(NEn,x,w,E1,e2)

    case ( CC_TYPE_G_LEGENDRE ) 

       call Gauss_Legendre_Rec(NEn, 0, E1, E2, x, w)

    case ( CC_TYPE_TANH_SINH ) 

       delta = 2.e-2_dp * abs(E2-E1) / real(NEn,dp)
       call TanhSinh_Exact(NEn, x, w, E1, E2, p=delta)

    case default
       call die('Could not determine the non-equilibrium line contour')
    end select

  end subroutine line_integral
  

  ! The residuals of the fermi-function at a real-energy
  subroutine nf_poles(PART,E, kT, Eta, Npol,c)

    use precision, only : dp
    use units, only : Pi

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: PART ! part of the contour
    real(dp), intent(in) :: E    ! at energy
    real(dp), intent(in) :: kT   ! temperature in Ry
    real(dp), intent(in) :: Eta  ! The lift of the equilibrium contour
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

    do while ( aimag(c(1)%c) < Eta ) 
       c(:)%c = c(:)%c + dcmplx(0._dp,2._dp*Pi*kT)
    end do

  end subroutine nf_poles
  

  

! ##################################################################
! ##      Generate (close to) real axis energy contour            ##
! ##           Transmission and DOS calculation                   ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
  subroutine transmission(E1,E2,nEq_Eta,NEn,contour)

    use precision, only : dp

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E1,E2    ! energy parameters 
    real(dp), intent(in) :: nEq_Eta    ! state broadening in Ry
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
       contour(ic)%c = dcmplx(E1+(ic-1)*delta, nEq_Eta)
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
  subroutine phonon(E1,E2,nEq_Eta,NEn,contour)

    use precision, only : dp

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E1,E2    ! energy parameters 
    real(dp), intent(in) :: nEq_Eta    ! state broadening in Ry
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
       contour(ic)%c = dcmplx(E1+(ic-1)*delta, nEq_Eta)**2
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
    c => contour(NEn_eq*2+1:NEn-N_transport)
  end function contour_nEq

  function contour_Transport() result(c)
    type(ts_ccontour), pointer :: c(:)
    c => contour(NEn-N_transport+1:NEn)
  end function contour_Transport

end module m_ts_contour
