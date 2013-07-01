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

  implicit none

  ! This module will also contain all the contour variables
  integer, save :: NEn_eq
  integer, save :: NEn  ! Number of energy points in the contour path
  integer, save :: PNEn ! Number of energy points in the contour path (divisible by Nodes)

  ! Contour path
  type(ts_ccontour), dimension(:), pointer, save :: contour
  type(ts_ccontour), dimension(:), pointer, save :: contourL => null()
  type(ts_ccontour), dimension(:), pointer, save :: contourR => null()
  type(ts_ccontour), dimension(:), pointer, save :: contour_neq => null()

  public :: NEn, PNEn, contour, contourL, contourR, contour_neq
  public :: setup_contour, io_contour, print_contour
  public :: sort_contour
  public :: init_Gauss_Fermi_plus_Line
  private

contains


! Routine for creating the contour
  subroutine setup_contour(IsVolt,C_eq_line,C_eq_circle,C_neq_tail, C_neq_mid, &
       EfL,Ef0,EfR, &
       NCircle,NLine,Npol,Nvolt, Nvolt_tail,NVolt_mid, &
       Emin,Emax,Ntransport, &
       CCEmin,GFEta,kT)
    
    use precision, only : dp
    use parallel,  only : IONode, Nodes, operator(.PARCOUNT.)
    use m_ts_aux, only : gaufermi0, gaufermi2

! **********************
! * INPUT variables    *
! **********************
    logical,  intent(in) :: IsVolt ! Do we have a volt
    ! The four governing methods
    integer, intent(in) :: C_eq_line, C_eq_circle, C_neq_tail, C_neq_mid
    real(dp), intent(in) :: EfL ! Left Fermi shift
    real(dp), intent(in) :: Ef0 ! equilibrium Fermi shift
    real(dp), intent(in) :: EfR ! Right Fermi shift
    ! The different contour path parts
    integer, intent(in)  :: Ncircle, Nline, Npol,Nvolt, Nvolt_tail, Nvolt_mid
    real(dp), intent(in) :: Emin ! Minimum transport energy
    real(dp), intent(in) :: Emax ! Maximum energy
    integer, intent(in)  :: Ntransport
    real(dp), intent(in) :: CCEmin, GFEta, kT

! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c(:) => null()
    integer :: NE_equilibrium
    integer :: i


    ! Determine how many contour points we have...  
    if ( IsVolt ) then
       NEn = 2*(Npol+Nline+Ncircle)+Nvolt+Ntransport
    else
       NEn = Npol+Nline+Ncircle+Ntransport
    end if
    PNEn = Nodes .PARCOUNT. NEn
    if ( NEn <= 0 ) then  ! if there are not any points we definetely have an error
       call die("ERROR: setup_contour: no energy points specified")
    end if
    
    ! Equilibrium contour points
    NE_equilibrium = Npol+Nline+Ncircle
    NEn_eq = NE_equilibrium

    ! Allocate the contour points
    nullify(contour)
    allocate(contour(NEn))
  
    ! Create the equilibrium contour in case we do not have 
    ! a voltage
    if ( (.not. IsVolt) .and. NE_equilibrium > 0 ) then

       c => contour(1:NE_equilibrium) 
       contourL => c
       call mod_HansSkriver(CC_PART_EQUI, &
            C_eq_line, C_eq_circle, &
            CCEmin, Ef0, kT, &
            Ncircle,Nline,Npol, c)
       
       ! Note we put a minus here because the integral we want is the
       ! negative of the pole-sum and C+L integral!!
       do i = 1 , NE_equilibrium
          c(i)%w = -c(i)%w
       end do

    else if ( NE_equilibrium > 0 ) then ! Do a voltage contour

       ! Do the left contour line
       c => contour(1:NE_equilibrium) 
       contourL => c
       call mod_HansSkriver(CC_PART_LEFT_EQUI, &
            C_eq_line, C_eq_circle, &
            CCEmin + EfL - Ef0, EfL, kT, &
            Ncircle,Nline,Npol, c)
       
       ! Note we put a minus here because the integral we want is the
       ! negative of the pole-sum and C+L integral!!
       do i = 1 , NE_equilibrium
          c(i)%w = -c(i)%w
       end do
       
       ! Do the right contour line
       c => contour(NE_equilibrium+1:2*NE_equilibrium) 
       contourR => c
       call mod_HansSkriver(CC_PART_RIGHT_EQUI, &
            C_eq_line, C_eq_circle, &
            CCEmin + EfR - Ef0, EfR, kT, &
            Ncircle,Nline,Npol, c)

       ! Note we put a minus here because the integral we want is the
       ! negative of the pole-sum and C+L integral!!
       do i = 1 , NE_equilibrium
          c(i)%w = -c(i)%w
       end do
       
       ! The voltage contour
       c => contour(2*NE_equilibrium+1:2*NE_equilibrium+NVolt) 
       contour_neq => c
       select case ( C_neq_tail ) ! check the bias contour
       case ( CC_TYPE_NEQ_SOMMERFELD ) ! 1. order
          
          call sommerfeld(CC_PART_NON_EQUI,EfR,EfL, &
               kT,GFeta, &
               NVolt, c)
          
       case ( CC_TYPE_NEQ_TAIL_G_NF_0kT )
          
          call Gauss_Fermi_plus_line(CC_PART_NON_EQUI,C_neq_tail,C_neq_mid, &
               gaufermi0, &
               EfR, EfL, 0._dp, &
               kT, GFeta, &
               NVolt, c, NVolt_tail, Nvolt_mid)

       case ( CC_TYPE_NEQ_TAIL_G_NF_2kT )
          
          call Gauss_Fermi_plus_line(CC_PART_NON_EQUI,C_neq_tail,C_neq_mid, &
               gaufermi2, &
               EfR, EfL, 2._dp*kT, &
               kT, GFeta, &
               NVolt, c, NVolt_tail, Nvolt_mid)

       case ( CC_TYPE_NEQ_TAIL_G_LAGUERRE )
          
          call die('not fully implemented')

       case ( CC_TYPE_NEQ_G_HERMITE )
          
          call die('not fully implemented')
          
       case default
          if ( IONode ) &
               write(*,*) 'ERROR: Contour not defined'
          call die('ERROR:  setup_contour: Contour not defined') 
       end select

    end if
    
    ! Finally we add the transport energy points
    if ( Ntransport > 0 ) then
       c => contour(NEn-Ntransport+1:NEn) 
       select case ( C_eq_line ) 
       case ( CC_TYPE_TRANSPORT ) 
          call transmission(Emin,Emax,GFeta,Ntransport, c)
       case ( CC_TYPE_TRANS_PHONON ) ! Phonon energy-points
          if (NEn /= Ntransport) then
             call die('ERROR: when doing phonon transport, only &
                  &use transport energy-points on the contour')
          end if
          call phonon(Emin,Emax,GFeta,Ntransport,c)
       case default
          call die('Unrecognized transport contour')
       end select
    end if

  end subroutine setup_contour




! ###################################################
! #                                                 #
! #       Different contour path variables          #
! # We define here the routines for generating      #
! # the different contour points.                   #
! #                                                 #
! ###################################################

!------------------------------------------------------------
!
!     Modified Hans Skriver Contour
!
  subroutine mod_HansSkriver(PART, C_eq_line, C_eq_circle, &
       E1,E2, kT, &
       Ncircle,Nline,Npol,contour)

    use precision, only : dp
    use parallel, only : IONode
    use units, only : Pi
    use m_ts_aux, only : nf
    use m_ts_aux, only : gaufermi10, gaufermi20
    use m_ts_cctype
    use m_gauss_quad, only : Gauss_Legendre_Rec, Gauss_Chebyshev_Exact

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: PART          ! part of the contour
    ! The methods used for the equilibrium contour
    integer,  intent(in) :: C_eq_line, C_eq_circle
    real(dp), intent(in) :: E1, E2        ! energy parameters 
    real(dp), intent(in) :: kT            ! temperature in Ry
    integer,  intent(in) :: Ncircle,Nline,Npol

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), target, intent(out) :: contour(Ncircle+Nline+Npol)

    ! Modified Hans Skriver:
    integer, parameter :: NT = 10   ! start line in modified HS at E2-NT*kT

! ***********************
! * LOCAL variables     *
! ***********************
    ! For temporary shift in the array (then we do not need an ic counter)
    type(ts_ccontour), pointer           :: c(:) => null()
    real(dp), dimension(:), allocatable :: theta,x,wt

    ! Various constants for calculating the contour
    real(dp) :: D, Delta, gamma
    real(dp) :: alpha, R, beta, y
    complex(dp) :: ztmp, z0

    ! loop variables
    integer :: i,j
    
    ! Parameters
    D     = E2-E1
    Delta = Npol*2.0d0*Pi*kT
    gamma = NT*kT
    alpha = dATAN(Delta/(D-gamma))
    R     = dsqrt(Delta*Delta + (D - gamma)*(D - gamma))/ &
         (2d0*cos(alpha))
    z0    = dcmplx(E1 + R, 0d0)
    beta  = dasin(Delta/R)

    ! Residuals
    c => contour(1:Npol)
    call nf_poles(PART, E2,kT,Npol,c)

    ! Line contour
    c => contour(Npol+1:Npol+Nline)
    allocate(wt(Nline),x(Nline))
    call memory('A','D',2*Nline,'mkCplxContour') 
    if(NT.eq.10) then
       call gaufermi10(Nline,x,wt)
    elseif(NT.eq.20) then
       call gaufermi20(Nline,x,wt)
    else
       if(IONode) then
          write(*,*) 'ERROR: ' 
          write(*,*) 'No Gauss quadrature for Fermi function '
       endif
       call die('No Gauss quadrature for Fermi function ')
    end if
    do i = 1 , Nline
       j = Nline-i+1         !reverse
       c(i)%c = dcmplx(x(j)*kT + E2,Delta)
       c(i)%w = -wt(j)*kT*dcmplx(1._dp,0._dp)
       c(i)%part = PART
       c(i)%type = C_eq_line
    end do
    call memory('D','D',2*Nline,'mkCplxContour')         
    deallocate(wt,x)

    ! Circle contour
    c => contour(Npol+Nline+1:Npol+Nline+Ncircle)
    allocate(wt(Ncircle),theta(Ncircle))
    call memory('A','D',2*Ncircle,'mkCplxContour')

    if ( C_eq_circle == CC_TYPE_EQ_CIRC_G_LEG ) then
       call Gauss_Legendre_Rec(Ncircle, 0, beta, Pi, theta, wt)
    else if ( C_eq_circle == CC_TYPE_EQ_CIRC_G_CH_O ) then
       ! We switch the integration bounds and change the sign of the weight
       ! this will give the correct path when out-put
       call Gauss_Chebyshev_Exact(Ncircle, theta, wt, a=Pi, b=beta, &
            method=1, pure=.false.)
       wt(:) = -wt(:)
    else if ( C_eq_circle == CC_TYPE_EQ_CIRC_G_CH_C ) then
       ! We switch the integration bounds and change the sign of the weight
       ! this will give the correct path when out-put
       call Gauss_Chebyshev_Exact(Ncircle, theta, wt, a=Pi, b=beta, &
            method=0, pure=.false.)
       wt(:) = -wt(:)
    end if

    do i = 1 , Ncircle 
       ztmp = R*exp(dcmplx(0._dp,theta(i)))
       c(i)%c = z0 + ztmp
       ! Factor i, comes from Ed\theta=dE=iR e^{i\theta}
       c(i)%w = wt(i)*nf((c(i)%c-E2)/kT)*dcmplx(0._dp,1._dp)*ztmp
       c(i)%part = PART
       c(i)%type = C_eq_circle
    end do

    call memory('D','D',2*Ncircle,'mkCplxContour') 
    deallocate(wt,theta)

  end subroutine mod_HansSkriver


! The residuals of the fermi-function at a real-energy
  subroutine nf_poles(PART,E, kT, &
       Npol,c)

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
       c(i)%type = CC_TYPE_EQ_RES
    end do
    
  end subroutine nf_poles



!------------------------------------------------------------
!
!     1. order Sommerfeld expansion using 2kT increments
!      
! We assume that the function to be integrated varies slowly on the
! kT-scale         
  subroutine sommerfeld(part,E1,E2, &
       kT,GFeta, &
       NEn,contour)

    use precision, only : dp
    use parallel, only : IONode
    use units, only : Pi

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: part          ! part of the contour
    real(dp), intent(in) :: E1, E2        ! energy parameters 
    real(dp), intent(in) :: kT            ! temperature in Ry
    real(dp), intent(in) :: GFeta         ! state broadening in Ry
    integer,  intent(in) :: NEn

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), intent(out) :: contour(NEn)

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
    
    if (E1 .gt. E2) then
       EE2=E1
       EE1=E2
    else
       EE2=E2
       EE1=E1
    end if
       
    delta = (EE2 - EE1)/(NEn-3) 
       
    etaSF = kT
    
    rtmp = kT*kT*Pi*Pi/(12.d0*etaSF)
    
    ! Denote the part
    contour(:)%part = part
    ! The sommerfeld must be the non-equilibrium contour part
    contour(:)%type = CC_TYPE_NEQ_SOMMERFELD

    ! Assign the contour points and weights
    contour(1)%c     = dcmplx(EE1 - etaSF,GFeta)
    contour(1)%w     = dcmplx(0.25d0*delta + rtmp,0d0)
    
    contour(2)%c     = dcmplx(EE1 + etaSF,GFeta)
    contour(2)%w     = dcmplx(0.25d0*delta - rtmp,0d0)
    
    contour(NEn-1)%c = dcmplx(EE2 - etaSF,GFeta)
    contour(NEn-1)%w = dcmplx(0.25d0*delta - rtmp,0d0)
    
    contour(NEn)%c   = dcmplx(EE2 + etaSF,GFeta)
    contour(NEn)%w   = dcmplx(0.25d0*delta + rtmp,0d0)
    
    do i = 3 , NEn - 2
       contour(i)%c = dcmplx(delta*(i-2) + EE1,GFeta) 
       contour(i)%w = dcmplx(delta,0d0)
    end do              !i
     
    if ( E1 .gt. E2 ) then
       do i = 1 , NEn
          contour(i)%w = -contour(i)%w
       end do
    end if                 !E1>E2
    
  end subroutine sommerfeld

!------------------------------------------------------------
!
!     Gaussian quadrature, using Gauss-Fermi for the tails of the non-equilibrium
!      
! We assume that the function to be integrated varies slowly on the
! kT-scale         
  subroutine Gauss_Fermi_plus_Line(part,G_F_method,line_method, GaussFermi, &
       E1,E2, gamma, &
       kT,GFeta, &
       NEn,contour,Ntail,Nmid)

    use units, only : Pi
    use precision, only : dp

    interface
       subroutine GaussFermi(n,x,w)
         use precision, only : dp
         integer, intent(in) :: n 
         real(dp), intent(out) :: x(n),w(n)
       end subroutine GaussFermi
    end interface

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: part          ! part of the contour
    integer,  intent(in) :: G_F_method    ! The type of tail method
    integer,  intent(in) :: line_method   ! The type of middle line method
    real(dp), intent(in) :: E1, E2        ! energy parameters 
    real(dp), intent(in) :: gamma         ! energy parameters for the temperature overlap
    real(dp), intent(in) :: kT            ! temperature in Ry
    real(dp), intent(in) :: GFeta         ! state broadening in Ry
    integer,  intent(in) :: NEn, Ntail, Nmid


! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), intent(out) :: contour(NEn)

! ***********************
! * LOCAL variables     *
! ***********************
    real(dp) :: wlt(20), xlt(20) ! Used to obtain the weights etc.
    ! For holding the energies:
    real(dp) :: EE1,EE2
    
    ! Loop variables
    integer :: i,j

    if ( NEn /= 2*Ntail + Nmid ) call die('Wrong number of points &
         for the tail and middle')

    if ( E1 > E2 ) then
       EE2 = E1
       EE1 = E2
    else
       EE2 = E2
       EE1 = E1
    end if
    
    call GaussFermi(Ntail,xlt,wlt)

    do i = 1 , Ntail
       j = Ntail - i + 1 ! reverse
       contour(i)%c          = dcmplx(-xlt(j)*kT + EE1,GFeta)
       contour(i)%w          = wlt(j)*kT
       contour(i)%part       = part
       contour(i)%type       = G_F_method
       contour(NEn+1-i)%c    = dcmplx(xlt(j)*kT + EE2,GFeta)
       contour(NEn+1-i)%w    = wlt(j)*kT
       contour(NEn+1-i)%part = part
       contour(NEn+1-i)%type = G_F_method
    end do

    ! Populate the middle line with the method asked for
    call nEq_line_quad(PART,line_method,EE1,EE2,gamma, &
         GFEta,kT,Nmid,contour(Ntail+1))
    
    if ( E1 > E2 ) then
       do i = 1 , NEn
          contour(i)%w = -contour(i)%w
       end do
    end if                 !E1>E2
    
  end subroutine Gauss_Fermi_plus_Line

  ! Determine the number of tail-points and 
  subroutine init_Gauss_Fermi_plus_Line(E1, E2, &
       kT, GFeta, NEn, Ntail, Nmid, GF_N_kT)
    
    use units, only : Pi
    use precision, only : dp
    use parallel, only : IONode

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E1, E2  ! energy parameters 
    real(dp), intent(in) :: kT      ! temperature in Ry
    real(dp), intent(in) :: GFeta   ! state broadening in Ry
    integer,  intent(in) :: NEn

! ***********************
! * OUTPUT variables    *
! ***********************
    integer,  intent(out) :: Ntail, Nmid, GF_N_kT

! ***********************
! * LOCAL variables     *
! ***********************

! Modified REAL AXIS CONTOUR:
    integer, parameter :: NGauF = 8 ! number of points [-inf,E2+NT*kT]
    integer, parameter :: NTGau = 2
    integer :: NGau     ! This holds how many points on the Fermi line
    
    ! Loop variables
    integer :: i,j

    ! first determine how many points to use for the fermiline
    GF_N_kT = NTGau
    if ( abs(E2) + abs(E1) .le. 3._dp*NTGau*kT ) then
       GF_N_kT = 0
    endif
    
    NGau = NEn*(2._dp+GF_N_kT)*kT/dabs(E1-E2) + 1._dp
    if ( NGau .gt. NGauF + GF_N_kT) then
       NGau = NGauF + GF_N_kT
    end if
    if(NEn .lt. 2*NGau+8) then
       NGau = (NEn-8)/2
    endif

    ! We have now determined the number of Gaussian quadrature points on the
    ! Fermi line

    if ( NGau .le. 0 ) then
       if (IONode) then
          write(*,*) &
               'ERROR: No. points=',NEn, &
               ' not valid for real axis gaussian quadrature, min=', 10
       end if
       call die('ERROR: Contour: too few points for real axis int') 
    end if

    Ntail = NGau
    Nmid  = NEn - 2*NGau

  end subroutine init_Gauss_Fermi_plus_Line
  

  subroutine nEq_line_quad(PART,TYPE,E1,E2,E_broad,GFEta,kT,NEn,contour)

    use precision, only : dp
    use units, only : Pi
    use m_ts_aux, only : nf1
    use m_gauss_quad ! Just all, many routines can be used

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: PART, TYPE
    real(dp), intent(in) :: E1,E2    ! energy parameters
    real(dp), intent(in) :: E_broad  ! This is the broadening from the fermi-tail.
!                                    ! Hence it is the integration part which is cut off from the line
    real(dp), intent(in) :: GFeta    ! state broadening in Ry
    real(dp), intent(in) :: kT       ! temperature in Ry
    integer,  intent(in) :: NEn      ! No. contour points

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), intent(out) :: contour(NEn)   ! points for contour

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: i, Ni
    real(dp) :: rtmp, delta
    real(dp), allocatable :: x(:), w(:)
    
    ! sort the energies
    if ( E1 > E2 ) call die('nEq-line: Input MUST be sorted for the energies')
    
    do i = 1 , NEn
       contour(i)%part = PART
       contour(i)%type = TYPE
    end do

    allocate(x(NEn),w(NEn))

    select case ( TYPE )
    case ( CC_TYPE_NEQ_MID_SIMP_EXT )
       ! in simpson we count: i = 0, ..., N
       Ni = NEn - 1
       
       delta = (E2 - E1 - 2._dp*E_broad)/real(Ni,dp)
       do i = 1 , NEn
          rtmp = E1 + E_broad + delta * (i-1)
          contour(i)%c = dcmplx(rtmp,GFeta) 
          contour(i)%w = dcmplx(delta,0._dp)* &
               (nf1((rtmp-E2)/kT) - nf1((rtmp-E1)/kT))
       end do

       ! extended simpsons rule
       contour(1    )%w = contour(1    )%w*17.d0/48.d0
       contour(2    )%w = contour(2    )%w*59.d0/48.d0
       contour(3    )%w = contour(3    )%w*43.d0/48.d0
       contour(4    )%w = contour(4    )%w*49.d0/48.d0
       contour(NEn-3)%w = contour(NEn-3)%w*49.d0/48.d0
       contour(NEn-2)%w = contour(NEn-2)%w*43.d0/48.d0
       contour(NEn-1)%w = contour(NEn-1)%w*59.d0/48.d0
       contour(NEn  )%w = contour(NEn  )%w*17.d0/48.d0

    case ( CC_TYPE_NEQ_MID_SIMP_COMP )
       ! in simpson we count: i = 0, ..., N
       Ni = NEn - 1
       if ( mod(Ni,2) /= 0 ) call die('Composite Simpson rule &
            &is only allowed for even N, please increase nEq points.')
       
       delta = (E2 - E1 - 2._dp*E_broad)/real(Ni,dp)
       do i = 1 , NEn
          rtmp = E1 + E_broad + delta * (i-1)
          contour(i)%c = dcmplx(rtmp,GFeta) 
          contour(i)%w = dcmplx(delta,0._dp)* &
               (nf1((rtmp-E2)/kT) - nf1((rtmp-E1)/kT)) / 3._dp
       end do
       
       ! Correct the weights for the composite simpson rule
       do i = 2 , NEn - 1, 2
          contour(i)%w = contour(i)%w * 4._dp
       end do
       do i = 3 , NEn - 2, 2
          contour(i)%w = contour(i)%w * 2._dp
       end do

    case ( CC_TYPE_NEQ_MID_SIMP_38 )
       ! in simpson we count: i = 0, ..., N
       Ni = NEn - 1
       
       delta = (E2 - E1 - 2._dp*E_broad)/real(Ni,dp)
       do i = 1 , NEn
          rtmp = E1 + E_broad + delta * (i-1)
          contour(i)%c = dcmplx(rtmp,GFeta) 
          contour(i)%w = dcmplx(delta,0._dp)* &
               (nf1((rtmp-E2)/kT) - nf1((rtmp-E1)/kT)) * 3._dp / 8._dp
       end do
       
       ! Correct the weights for the 3/8 Simpson
       call die('currently not fully implemented')
    case ( CC_TYPE_NEQ_MID_MID )
       ! set boundaries for gaussian quadrature
       delta = (E2 - E1 - 2._dp*E_broad)/real(NEn,dp)
       do i = 1 , NEn
          ! move into the middle of the current segment
          rtmp = E1 + E_broad + delta * ( real(i,dp)-.5_dp )
          contour(i)%c = dcmplx(rtmp,GFeta) 
          contour(i)%w = dcmplx(delta,0._dp)* &
               (nf1((rtmp-E2)/kT) - nf1((rtmp-E1)/kT))
       end do

    case ( CC_TYPE_NEQ_MID_G_LEG ) 
       call Gauss_Legendre_Rec(NEn, 0, E1+E_broad, E2-E_broad, x, w)
       do i = 1 , NEn
          contour(i)%c = dcmplx(x(i),GFeta)
          contour(i)%w = dcmplx(w(i),0._dp)
       end do

    case ( CC_TYPE_NEQ_MID_G_CH_O ) 
       ! interchange bounds to correct out-put
       call Gauss_Chebyshev_Exact(NEn, x, w, &
            b=E1+E_broad, a=E2-E_broad, method=1, pure=.false.)
       do i = 1 , NEn
          contour(i)%c = dcmplx(x(i),GFeta)
          contour(i)%w = dcmplx(-w(i),0._dp)
       end do

    case ( CC_TYPE_NEQ_MID_G_CH_C ) 
       ! interchange bounds to correct out-put
       call Gauss_Chebyshev_Exact(NEn, x, w, &
            b=E1+E_broad, a=E2-E_broad, method=0, pure=.false.)
       do i = 1 , NEn
          contour(i)%c = dcmplx(x(i),GFeta)
          contour(i)%w = dcmplx(-w(i),0._dp)
       end do

    case default
       call die('Could not determine the non-equilibrium line contour')
    end select

    deallocate(x,w)

  end subroutine nEq_line_quad

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

end module m_ts_contour
