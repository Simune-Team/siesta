!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_contour
!
! Routines that are used to setup the contour for integration of the GFs
! 
! Note: I am unsure whether to move the METHOD variables to the m_ts_cctype
!       as they "belong" there. However, the relavant routines are located
!       here.
!
! A couple of routines are included in this module
!   1) setup_contour
!   2) io_contour
!   3) print_contour

! Furthermore we have a couple of local routines 
! used for generating the contour points
!   4) mod_HansSkriver
!   5) sommerfeld
!   6) gaussian_quadrature


! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_cctype

  implicit none

! The different available contours
  integer, parameter, public :: CC_METHOD_MOD_HANSSKRIVER = 1
  integer, parameter, public :: CC_METHOD_SOMMERFELD      = 2
  integer, parameter, public :: CC_METHOD_GAUSSFERMI      = 3

! This module will also contain all the contour variables
  integer, save :: NEn  ! Number of energy points in the contour path
  integer, save :: PNEn ! Number of energy points in the contour path (divisible by Nodes)

! Contour path
  type(ts_ccontour), dimension(:), pointer, save :: contour

  public :: NEn, PNEn, contour
  public :: setup_contour, io_contour, print_contour

  private

contains


! Routine for creating the contour
  subroutine setup_contour(IsVolt,Cmethod,EfL,Ef0,EfR, &
       NCircle,NLine,Npol,Nvolt,Emin,Emax,Ntransport, &
       CCEmin,GFEta,kT)
    
    use precision, only : dp
    use sys,       only : die
    use parallel,  only : IONode, Nodes, operator(.PARCOUNT.)

! **********************
! * INPUT variables    *
! **********************
    logical,  intent(in) :: IsVolt ! Do we have a volt
    integer, intent(in)  :: cmethod ! The method of the voltage contour
    real(dp), intent(in) :: EfL ! Left Fermi shift
    real(dp), intent(in) :: Ef0 ! equilibrium Fermi shift
    real(dp), intent(in) :: EfR ! Right Fermi shift
! The different contour path parts
    integer, intent(in)  :: Ncircle, Nline, Npol,Nvolt
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

    ! Allocate the contour points
    nullify(contour)
    allocate(contour(NEn))
  
! Create the equilibrium contour in case we do not have 
! a voltage
    if ( .not. IsVolt .and. NE_equilibrium > 0 ) then

       c => contour(1:NE_equilibrium) 
       call mod_HansSkriver(CC_PART_EQUI, &
            CCEmin, Ef0, &
            kT,GFeta, &
            Ncircle,Nline,Npol, c)
              
! Note we put a minus here because the integral we want is the
! negative of the pole-sum and C+L integral!!
       do i = 1 , NE_equilibrium
          c(i)%w = -c(i)%w
       end do

    else if ( NE_equilibrium > 0 ) then ! Do a voltage contour

! Do the left contour line
       c => contour(1:NE_equilibrium) 
       call mod_HansSkriver(CC_PART_LEFT_EQUI, &
            CCEmin + EfL - Ef0, EfL, &
            kT,GFeta, &
            Ncircle,Nline,Npol, c)
       
! Note we put a minus here because the integral we want is the
! negative of the pole-sum and C+L integral!!
       do i = 1 , NE_equilibrium
          c(i)%w = -c(i)%w
       end do
       
! Do the right contour line
       c => contour(NE_equilibrium+1:2*NE_equilibrium) 
       call mod_HansSkriver(CC_PART_RIGHT_EQUI, &
            CCEmin + EfR - Ef0, EfR, &
            kT,GFeta, &
            Ncircle,Nline,Npol, c)

! Note we put a minus here because the integral we want is the
! negative of the pole-sum and C+L integral!!
       do i = 1 , NE_equilibrium
          c(i)%w = -c(i)%w
       end do
       
! The voltage contour
       c => contour(2*NE_equilibrium+1:2*NE_equilibrium+NVolt) 
       if (Cmethod.eq.CC_METHOD_SOMMERFELD) then ! 1. order
          
          call sommerfeld(CC_PART_NON_EQUI,EfR,EfL, &
               kT,GFeta, &
               NVolt, c)
          
       else if(Cmethod.eq.CC_METHOD_GAUSSFERMI) then
          
          call gaussian_quadrature(CC_PART_NON_EQUI,EfR,EfL, &
               kT,GFeta, &
               NVolt, c)

       else
          if ( IONode ) &
               write(*,*) 'ERROR: Contour not defined'
          call die('ERROR:  setup_contour: Contour not defined') 
       end if

    end if
    
! Finally we add the transport energy points
    if ( Ntransport > 0 ) then
       c => contour(NEn-Ntransport+1:NEn) 
       call transmission(Emin,Emax, &
            GFeta,Ntransport, c)
    end if

  end subroutine setup_contour


  
  ! Routine for "pretty" printing the contour points out in the out file
  subroutine print_contour()
    use parallel, only : IONode
    use units,    only : eV

! **********************
! * LOCAL variables    *
! **********************
    type(ts_ccontour), pointer :: c
    character(len=6) :: ctype
    integer :: i, part, type
    
    ! Initialize variables
    part = -1
    type = -1
    nullify(c)

    if ( IONode ) then
       write(*,'(a)') "transiesta: contour integration path:"
       write(*,'(1x,a6,''   '',2(tr1,a12),2(tr1,a14))') &
            "Type  ","Re(c)[eV]","Im(c)[eV]","Re(weight)","Im(weight)"
       do i = 1 , NEn
          ! loop !
          c => contour(i)
          if ( part /= c%part ) then
             if ( c%part == CC_PART_EQUI ) then
                write(*,'(1x,a)') "Equilibrium:"
             else if ( c%part == CC_PART_LEFT_EQUI ) then
                write(*,'(1x,a)') "Left equilibrium:"
             else if ( c%part == CC_PART_RIGHT_EQUI ) then
                write(*,'(1x,a)') "Right equilibrium:"
             else if ( c%part == CC_PART_NON_EQUI ) then
                write(*,'(1x,a)') "Non-equilibrium:"
             else if ( c%part == CC_PART_TRANSPORT ) then
                write(*,'(1x,a)') "Transport:"
             end if
             part = c%part
          end if
          if ( type /= c%type ) then
             if ( c%type == CC_TYPE_RES ) then
                ctype = 'resi'
             else if ( c%type == CC_TYPE_FERMI ) then
                ctype = 'fermi'
             else if ( c%type == CC_TYPE_CIRCLE ) then
                ctype = 'circle'
             else if ( c%type == CC_TYPE_NON_EQUI ) then
                ctype = 'noneq'
             else if ( c%type == CC_TYPE_GAUSS_FERMI ) then
                ctype = 'gaussF'
             else if ( c%type == CC_TYPE_GAUSS_QUAD ) then
                ctype = 'gaussQ'
             else if ( c%type == CC_TYPE_TRANSPORT ) then
                ctype = 'trans'
             end if
             type = c%type
          end if
          
          ! Write out the contour information:
          write(*,'(1x,a6,'' : '',tr1,2(f12.5,tr1),2(f14.9,tr1))') &
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
    character(len=len_trim(slabel)+5) :: fname
    integer :: i, unit, part

    if ( IONode ) then
       fname = trim(slabel) // '.TSCC'
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
  subroutine mod_HansSkriver(part,E1,E2, &
       kT,GFeta, &
       Ncircle,Nline,Npol,contour)

    use precision, only : dp
    use parallel, only : IONode
    use sys, only : die
    use units, only : Pi
    use m_ts_aux_rout, only : nf
    use m_ts_aux_rout, only : gaufermi10, gaufermi20, gauss

! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: part          ! part of the contour
    real(dp), intent(in) :: E1, E2        ! energy parameters 
    real(dp), intent(in) :: kT            ! temperature in Ry
    real(dp), intent(in) :: GFeta         ! state broadening in Ry
    integer,  intent(in) :: Ncircle,Nline,Npol

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_ccontour), target, intent(out) :: contour(Ncircle+Nline+Npol)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
    real(dp) :: alpha, R, beta
    complex(dp) :: ztmp, z0

! loop variables
    integer :: i,j

!               
!     Parameters
!     
    D     = E2-E1
    Delta = Npol*2.0d0*Pi*kT
    gamma = NT*kT
    alpha = dATAN(Delta/(D-gamma))
    R     = dsqrt(Delta*Delta + (D - gamma)*(D - gamma))/ &
         (2d0*cos(alpha))
    z0    = dcmplx(E1 + R, 0d0)
    beta  = dasin(Delta/R)

!
!     Residuals:
!        
    c => contour(1:Npol)
    do i = 1 , Npol 
       c(i)%c = dcmplx(E2,Pi*kT*(2.0d0*(i-1)+1d0))
       c(i)%w = dcmplx(0d0,2d0*Pi*kT) 
       c(i)%part = part
       c(i)%type = CC_TYPE_RES
    end do                 !i

!
!     Line contour:
!        
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
          write(*,*) &
               'No Gauss quadrature for Fermi function '
       endif
       call die('No Gauss quadrature for Fermi function ')
    end if
    do i = 1 , Nline
       j = Nline-i+1         !reverse
       c(i)%c = dcmplx(x(j)*kT + E2,Delta)
       c(i)%w = -wt(j)*kT*dcmplx(1d0,0d0)
       c(i)%part = part
       c(i)%type = CC_TYPE_FERMI
    end do
    call memory('D','D',2*Nline,'mkCplxContour')         
    deallocate(wt,x)

!
!     Circle contour:
!
    c => contour(Npol+Nline+1:Npol+Nline+Ncircle)
    allocate(wt(Ncircle),theta(Ncircle))
    call memory('A','D',2*Ncircle,'mkCplxContour') 
    call gauss(Ncircle, 0, beta, Pi, theta, wt)       

    do i = 1 , Ncircle 
       ztmp = R*exp(theta(i)*dcmplx(0d0,1d0))
       c(i)%c = z0 + ztmp
       c(i)%w = wt(i)*nf((c(i)%c-E2)/(kT))*dcmplx(0d0,1d0)*ztmp
       c(i)%part = part
       c(i)%type = CC_TYPE_CIRCLE
    end do                 !i
    call memory('D','D',2*Ncircle,'mkCplxContour') 
    deallocate(wt,theta)

  end subroutine mod_HansSkriver



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
    use sys, only : die
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
    contour(:)%type = CC_TYPE_NON_EQUI

    ! Assign the contour points and weights
    contour(1)%c    = dcmplx(EE1 - etaSF,GFeta)
    contour(1)%w    = dcmplx(0.25d0*delta + rtmp,0d0)
    
    contour(2)%c    = dcmplx(EE1 + etaSF,GFeta)
    contour(2)%w    = dcmplx(0.25d0*delta - rtmp,0d0)
    
    contour(NEn-1)%c = dcmplx(EE2 - etaSF,GFeta)
    contour(NEn-1)%w = dcmplx(0.25d0*delta - rtmp,0d0)
    
    contour(NEn)%c = dcmplx(EE2 + etaSF,GFeta)
    contour(NEn)%w = dcmplx(0.25d0*delta + rtmp,0d0)
    
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
!     Gaussian quadrature, using gaufermi in the ends
!      
! We assume that the function to be integrated varies slowly on the
! kT-scale         
  subroutine gaussian_quadrature(part,E1,E2, &
       kT,GFeta, &
       NEn,contour)

    use parallel, only : IONode
    use sys, only : die
    use units, only : Pi
    use precision, only : dp
    use m_ts_aux_rout, only : nf1, gaufermi0, gaufermi2

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

! Modified REAL AXIS CONTOUR:
    integer, parameter :: NGauF = 8 ! number of points [-inf,E2+NT*kT]
    integer, parameter :: NTGau = 2
    real(dp) :: wlt(NGauF),xlt(NGauF) ! Used to obtain the weights etc.
    integer :: NGau ! This holds how many points on the Fermi line
    integer :: NTGauUse ! Intermediate for determining the Fermi line points
! For holding the energies:
    real(dp) :: EE1,EE2
! Different uses
    real(dp) :: rtmp, gamma, delta
    
! Loop variables
    integer :: i,j, Ni

    ! first determine how many points to use for the fermiline
    NTGauUse = NTGau
    if ( dabs(E2 - E1) .le. 3.d0*NTGau*kT ) then
       NTGauUse = 0
    endif
    
    NGau = NEn*(2.d0+NTGauUse)*kT/dabs(E1-E2) +1.0d0
    if ( NGau .gt. NGauF+NTGauUse) then
       NGau = NGauF + NTGauUse
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
               ' not valid for real axis gaussian quadrature, min', 10
       end if
       call die('ERROR: Contour: too few points for real axis int') 
    end if

    if ( E1 .gt. E2 ) then
       EE2 = E1
       EE1 = E2
    else
       EE2 = E2
       EE1 = E1
    end if
    
    if ( NTGauUse .eq. 2 ) then
       call gaufermi2(NGau,xlt,wlt)
    else if ( NTGauUse .eq. 0 ) then
       call gaufermi0(NGau,xlt,wlt)
    else
       if ( IONode ) then
          write(*,*) 'ERROR: ' 
          write(*,*) &
               'No Gauss quadrature for Fermi function '
       endif
       call die ('No Gauss quadrature for Fermi function')
    end if

    do i = 1 , NGau
       j = NGau - i + 1 ! reverse
       contour(i)%c       = dcmplx(-xlt(j)*kT + EE1,GFeta)
       contour(i)%w       = wlt(j)*kT*dcmplx(1d0,0d0)
       contour(i)%part    = part
       contour(i)%type    = CC_TYPE_GAUSS_FERMI
       contour(NEn+1-i)%c = dcmplx(xlt(j)*kT + EE2,GFeta)
       contour(NEn+1-i)%w = wlt(j)*kT*dcmplx(1d0,0d0)
       contour(NEn+1-i)%part = part
       contour(NEn+1-i)%type = CC_TYPE_GAUSS_FERMI
    end do

    gamma = NTGauUse*kT

! set boundaries for gaussian quadrature
    delta = (EE2 - EE1-2.*gamma)/(NEn-2*NGau-1)
    do i = NGau+1,NEn-NGau
       rtmp = delta*(i-NGau-1) + EE1+gamma
       contour(i)%c = dcmplx(rtmp,GFeta) 
       contour(i)%w = dcmplx(delta,0d0)* &
            (nf1((rtmp-EE2)/kT) -nf1((rtmp-EE1)/kT))
       contour(i)%part = part
       contour(i)%type = CC_TYPE_GAUSS_QUAD
    end do              !i

! extended simpsons rule
    contour(NGau+1)%w = contour(NGau+1)%w*17.d0/48.d0
    contour(NGau+2)%w = contour(NGau+2)%w*59.d0/48.d0
    contour(NGau+3)%w = contour(NGau+3)%w*43.d0/48.d0
    contour(NGau+4)%w = contour(NGau+4)%w*49.d0/48.d0
    Ni = NEn-NGau
    contour(Ni+1-1)%w = contour(Ni+1-1)%w*17.d0/48.d0
    contour(Ni+1-2)%w = contour(Ni+1-2)%w*59.d0/48.d0
    contour(Ni+1-3)%w = contour(Ni+1-3)%w*43.d0/48.d0
    contour(Ni+1-4)%w = contour(Ni+1-4)%w*49.d0/48.d0
 
    if ( E1 .gt. E2 ) then
       do i = 1 , NEn
          contour(i)%w = -contour(i)%w
       end do
    end if                 !E1>E2

  end subroutine gaussian_quadrature


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
    delta = (E2-E1)/real(NEn, dp)

    do ic = 1 , NEn
       contour(ic)%c = dcmplx(E1+(ic-0.5_dp)*delta, GFeta)
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


end module m_ts_contour
