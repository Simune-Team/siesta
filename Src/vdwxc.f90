! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2008.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.

MODULE m_vdwxc

! Implements the van der Waals functional kernel of
! M.Dion et al, PRL 92, 246401 (2004)
! Written by J.M.Soler. July 2007.

  use precision,   only: dp
  use mesh1D,      only: get_mesh, get_n, set_mesh, set_interpolation, &
                         derivative, integral
  use m_radfft,    only: radfft
  use flib_spline, only: generate_spline
  use m_recipes,   only: spline, splint

! BEGIN DEBUG
  use m_debug,     only: udebug     ! File unit for debug output
!  use plot_module, only: plot
! END DEBUG

  implicit none

! Called by xc routines
PUBLIC :: vdw_decusp, vdw_theta, vdw_get_qmesh, vdw_phi, vdw_set_kcut

! Called by debugging programs
 PUBLIC :: phiofr, phi_interp, phi_soft, phi_val, pofq, qofrho

PRIVATE  ! Nothing is declared public beyond this point

!  integer, parameter:: dp = kind(1.d0)

  ! Precision parameters for the integral defining phi in routine phi_val
  real(dp),parameter:: acutmin = 10.0_dp   ! Upper integration limit
  real(dp),parameter:: acutbyd  = 30.0_dp  ! Upper integration limit / d
  real(dp),parameter:: damax = 0.5_dp      ! Max. integration interval
  real(dp),parameter:: damin = 1.e-2_dp    ! Min. integration interval
  real(dp),parameter:: dabyd = 0.1_dp      ! Min. integration interval / d

  ! Precision parameter for the integral in routine dphi
  real(dp),parameter:: ashortbyd = 2.5_dp  ! Shorter integration limit / d

  ! Parameters for phi_soft
  ! (dsoft,phi0_soft)=(0.5,0.8)|(0.7|0.6)|(1.0|0.4)|(1.5,0.22)|(2.0,0.12)
  real(dp),parameter:: dsoft = 1.0_dp       ! Softening matching radius
  real(dp),parameter:: phi0_soft = 0.40_dp  ! phi_soft(0,0) (depends on dsoft)
  real(dp),parameter:: dd = 0.01_dp         ! Delta(d) for derivatives

  ! Mesh parameters for table phi(d1,d2)
  integer, parameter:: nd = 20               ! Number of d mesh points
  real(dp),parameter:: dcut = 30.0_dp        ! Max. value of d mesh
  real(dp),parameter:: ddmaxddmin = 20.0_dp  ! Last d mesh interval / first one

  ! Use routine dphi for better efficiency in setting phi table?
  logical,parameter:: use_dphi =.true. !(.true.=>efficiency|.false.=>accuracy)

  ! Set derivation methods to use for interpolation table
  character(len=*),parameter:: deriv_method = 'numeric'  !('numeric'|'interp')
  character(len=*),parameter:: interp_method= 'Spline' !('Lagrange'|'Spline')

  ! Parameters to find numerical derivatives of phi, to set interpolation
  real(dp),parameter:: ddbyd = 0.01_dp       ! Delta to find phi derivs / d
  real(dp),parameter:: ddmin = 0.001_dp      ! Min. delta to find phi derivs

  ! Mesh parameters for table of phi(q1,q2,r) and its Fourier transform
  integer, parameter:: nr = 1024             ! Radial mesh points (power of 2)
  integer, parameter:: mq = 20               ! Total number of q mesh points
  integer, parameter:: nq = mq-1             ! Effective number of q mesh points
  real(dp),parameter:: qcut = 5.0_dp         ! Max. value of q mesh
  real(dp),parameter:: dqmaxdqmin = 20.0_dp  ! Last q mesh interval / first one
  real(dp),parameter:: rcut = 100._dp        ! Radial cutoff: r>rcut => phi=0
  real(dp),parameter:: rmin = 1.e-6_dp       ! Min. radius as denominator
  real(dp),parameter:: rsoft = 0.0_dp        ! Soften kernel in r<rsoft

  ! Parameters for cutoff function, used in radial Fourier transforms of phi
  integer, parameter:: ncut1 =  8      ! cutoff(x)=(1-x**ncut1)**ncut2
  integer, parameter:: ncut2 =  4

  ! Parameters for saturate function, used to enforce that q<qcut
  integer, parameter:: nsat  = 12      ! xsat(x)=1-exp(-sum_n=1:nsat x**n/n)

  ! Parameters for saturate_inverse function
  real(dp),parameter:: xmaxbyxc = 100._dp   ! qmax/qcut
  real(dp),parameter:: ytol = 1.e-15_dp     ! Tol. for saturated q

  ! Private module variables and arrays
  real(dp):: dmesh(nd)                ! Mesh points for phi(d1,d2) table
  real(dp):: qmesh(mq)                ! Mesh points for phi(q1,q2,r)
  real(dp):: phi_table(0:3,0:3,nd,nd) ! Coefs. for bicubic interpolation
  logical :: phi_table_set=.false.    ! Has phi_table been set?
  logical :: qmesh_set=.false.        ! Has qmesh been set?
  logical :: kcut_set=.false.         ! Has kcut been set?
  real(dp):: phir(0:nr,mq,mq)         ! Table of phi(r)
  real(dp):: d2phidr2(0:nr,mq,mq)     ! Table of d2_phi/dr2
  real(dp):: dr                       ! r-mesh interval
  real(dp):: phik(0:nr,mq,mq)         ! Table of phi(k)
  real(dp):: d2phidk2(0:nr,mq,mq)     ! Table of d2_phi/dk2
  real(dp):: dk                       ! k-mesh interval
  real(dp):: kcut                     ! Planewave cutoff: k>kcut => phi=0
  integer :: nk                       ! # k points within kcut

!  real(dp):: dqmaxdqmin, qcut

CONTAINS

! -----------------------------------------------------------------------------

SUBROUTINE bcucof( n1, n2, x1, x2, y, dydx1, dydx2, d2ydx1dx2, c )
! Finds coefficients for bicubic interpolation
! Adapted from Numerical recipes
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1, n2
  REAL(dp),INTENT(IN) :: x1(n1)
  REAL(dp),INTENT(IN) :: x2(n2)
  REAL(dp),INTENT(IN) :: y(n1,n2)
  REAL(dp),INTENT(IN) :: dydx1(n1,n2)
  REAL(dp),INTENT(IN) :: dydx2(n1,n2)
  REAL(dp),INTENT(IN) :: d2ydx1dx2(n1,n2)
  REAL(dp),INTENT(OUT):: c(0:3,0:3,n1,n2)

  INTEGER  :: i1, i11, i12, i13, i14, i2, i21, i22, i23, i24
  REAL(dp) :: dx1, dx2, wt(16,16), z(16)
  DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
    8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
    2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
    2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
    -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
    -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
    -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

! Set coefs. for i1<n1 and i2<n2
  do i2 = 1,n2-1
    do i1 = 1,n1-1
      dx1 = x1(i1+1) - x1(i1)
      dx2 = x2(i2+1) - x2(i2)
      i11 = i1
      i12 = i1+1
      i13 = i1+1
      i14 = i1
      i21 = i2
      i22 = i2
      i23 = i2+1
      i24 = i2+1
      z( 1) = y(i11,i21)
      z( 2) = y(i12,i22)
      z( 3) = y(i13,i23)
      z( 4) = y(i14,i24)
      z( 5) = dydx1(i11,i21) * dx1
      z( 6) = dydx1(i12,i22) * dx1
      z( 7) = dydx1(i13,i23) * dx1
      z( 8) = dydx1(i14,i24) * dx1
      z( 9) = dydx2(i11,i21) * dx2
      z(10) = dydx2(i12,i22) * dx2
      z(11) = dydx2(i13,i23) * dx2
      z(12) = dydx2(i14,i24) * dx2
      z(13) = d2ydx1dx2(i11,i21) * dx1 * dx2
      z(14) = d2ydx1dx2(i12,i22) * dx1 * dx2
      z(15) = d2ydx1dx2(i13,i23) * dx1 * dx2
      z(16) = d2ydx1dx2(i14,i24) * dx1 * dx2
      z = matmul(wt,z)
      c(0:3,0:3,i1,i2) = reshape(z,(/4,4/),order=(/2,1/))
    end do ! i1
  end do ! i2

! Set c for i1=n1 and i2=n2 (valid only at the border)
  c(:,:,n1,:) = 0
  c(:,:,:,n2) = 0
  c(0,0,n1,:) = y(n1,:)
  c(0,0,:,n2) = y(:,n2)

END SUBROUTINE bcucof

! -----------------------------------------------------------------------------

function cutoff( x )

  implicit none
  real(dp),intent(in):: x
  real(dp)           :: cutoff

  if (x<=0._dp) then
    cutoff = 1
  else if (x>=1._dp) then
    cutoff = 0
  else
    cutoff = (1-x**ncut1)**ncut2
  end if

end function cutoff

!-----------------------------------------------------------------------------

real(dp) function dphi( d1, d2 )

! Finds phi(d1,d2) - phi(dmax,dmax), with dmax=max(d1,d2)

  real(dp),intent(in) :: d1, d2

  integer  :: ia, ib, n, nmesh, nshort
  real(dp) :: a, acut, b, da1, dan, deff, dmax, dmin, gamma, &
              pi, t, t0, w
  real(dp),allocatable :: amesh(:), c(:), dphida(:), dphidb(:), &
                          s(:), nu0(:), nu1(:), nu2(:)

! Find integration mesh
  dmax = max( abs(d1), abs(d2) )
  dmin = min( abs(d1), abs(d2) )
  deff = dmax
!  deff = sqrt(d1**2+d2**2)
  acut = max( acutmin, acutbyd * deff )
  da1 = min( damax, dabyd * deff )
  da1 = max( da1, damin )
  dan = damax
  n = get_n( 0._dp, acut, da1, dan )
  allocate( amesh(n), c(n), dphida(n), dphidb(n), s(n), &
            nu0(n), nu1(n), nu2(n) )
  call set_mesh( n, xmax=acut, dxndx1=dan/da1 )
  call get_mesh( n, nmesh, amesh )

! Find limit of shorter mesh
  nshort = n
  do ia = n-1,1,-1
    if (amesh(ia) > ashortbyd*deff) nshort = ia
  end do

! Find cos(a), sin(a), nu1(a), and nu2(a)
  pi = acos(-1._dp)
  gamma = 4*pi/9
  do ia = 1,n
    a = amesh(ia)
    c(ia) = cos(a)
    s(ia) = sin(a)
    if (d1==0._dp) then
      nu1(ia) = a*a / 2
    else
      nu1(ia) = a*a / 2 / (1-exp(-gamma*(a/d1)**2))
    end if
    if (d2==0._dp) then
      nu2(ia) = a*a / 2
    else
      nu2(ia) = a*a / 2 / (1-exp(-gamma*(a/d2)**2))
    end if
    if (dmax<=0._dp) then
      nu0(ia) = a*a / 2
    else
      nu0(ia) = a*a / 2 / (1-exp(-gamma*(a/dmax)**2))
    end if
  end do

! Make integral on variable a
  dphida(1) = 0
  do ia = 2,nshort
    a = amesh(ia)

    ! Make integral on variable b
    dphidb(1) = 0
    do ib = 2,n
      b = amesh(ib)

      w = 2*( (3-a*a)*b*c(ib)*s(ia) + (3-b*b)*a*c(ia)*s(ib) &
            + (a*a+b*b-3)*s(ia)*s(ib) - 3*a*b*c(ia)*c(ib) )/(a*b)**3

      t = 0.5_dp * ( 1/(nu1(ia)+nu1(ib)) + 1/(nu2(ia)+nu2(ib)) ) &
                 * ( 1/(nu1(ia)+nu2(ia))/(nu1(ib)+nu2(ib)) &
                   + 1/(nu1(ia)+nu2(ib))/(nu2(ia)+nu1(ib)) )

      t0 = 0.5_dp * ( 1/(nu0(ia)+nu0(ib)) + 1/(nu0(ia)+nu0(ib)) ) &
                  * ( 1/(nu0(ia)+nu0(ia))/(nu0(ib)+nu0(ib)) &
                    + 1/(nu0(ia)+nu0(ib))/(nu0(ia)+nu0(ib)) )

      dphidb(ib) = a*a * b*b * w * (t-t0)

    end do ! ib
    dphida(ia) = integral( n, dphidb ) - integral( ia, dphidb )

  end do ! ia
  dphi = 2/pi**2 * 2*integral( nshort, dphida )

  deallocate( amesh, c, dphida, dphidb, s, nu0, nu1, nu2 )

! BEGIN DEBUG
!  print'(a,2f12.6,i8,f12.6)', 'dphi: d1,d2, na, phi =', d1, d2, n, dphi
! END DEBUG

end function dphi

! -----------------------------------------------------------------------------

function dphi_fast( d1, d2 )

! Finds phi(d1,d2)-phi(dmax,dmax), with dmax=max(d1,d2), by 
!  - Direct integration if d < dsoft, where d=sqrt(d1**2+d2**2)
!  - Interpolation of phi_table if d > dsoft

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: dphi_fast

  real(dp):: d, dmax

  if (.not.phi_table_set) call set_phi_table()

  d = sqrt( d1**2 + d2**2 )
  dmax = max( d1, d2 )

  if (d < dsoft) then
    dphi_fast = dphi( d1, d2 )
  else
    dphi_fast = phi_interp( d1, d2 ) - phi_interp( dmax, dmax )
  end if

end function dphi_fast

! -----------------------------------------------------------------------------

function dphi_soft( d1, d2 )

! Finds phi_soft(d1,d2)-phi_soft(dmax,dmax), with dmax=max(d1,d2)

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: dphi_soft

  real(dp):: d, dmax

  d = sqrt( d1**2 + d2**2 )
  dmax = max( d1, d2 )

  if (d < dsoft) then
    dphi_soft = phi_soft( d1, d2 ) - phi_soft( dmax, dmax )
  else
    dphi_soft = dphi( d1, d2 )
  end if

end function dphi_soft

!-----------------------------------------------------------------------------

integer function iofd( d )

! Finds index i such that dmesh(i) <= d < dmesh(i+1)

  implicit none
  real(dp), intent(in) :: d

  real(dp),parameter :: amin = 1.e-12_dp
  real(dp),save:: a, b
  logical, save:: first_call = .true.

  if (first_call) then
    a = log( (dmesh(nd)-dmesh(nd-1)) / (dmesh(2)-dmesh(1)) ) / (nd-2)
    a = max( a, amin )
    b = (dmesh(2) - dmesh(1)) / (exp(a) - 1)
    first_call = .false.
  end if

  iofd = 1 + log( 1 + (d-dmesh(1))/b ) / a
  iofd = max( 1, iofd )
  iofd = min( nd-1, iofd )

end function iofd

!-----------------------------------------------------------------------------

integer function iofq( q )

! Finds index i such that qmesh(i) <= q < qmesh(i+1)

  implicit none
  real(dp), intent(in) :: q

  real(dp),parameter :: amin = 1.e-12_dp
  real(dp),save:: a, b
  logical, save:: first_call = .true.

  if (first_call) then
    a = log( (qmesh(mq)-qmesh(mq-1)) / (qmesh(2)-qmesh(1)) ) / (mq-2)
    a = max( a, amin )
    b = (qmesh(2) - qmesh(1)) / (exp(a) - 1)
    first_call = .false.
  end if

  iofq = 1 + log( 1 + (q-qmesh(1))/b ) / a
  iofq = max( 1, iofq )
  iofq = min( mq-1, iofq )

end function iofq

! -----------------------------------------------------------------------------

function phi_fast( d1, d2 )

! Finds hard phi(d1,d2) kernel by 
!  - Direct integration if d < dsoft, where d=sqrt(d1**2+d2**2)
!  - Interpolation of phi_table if d > dsoft

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: phi_fast

  real(dp):: d

  if (.not.phi_table_set) call set_phi_table()

  d = sqrt( d1**2 + d2**2 )

  if (d < dsoft) then
    phi_fast = phi_val( d1, d2 )
  else
    phi_fast = phi_interp( d1, d2 )
  end if

end function phi_fast

! -----------------------------------------------------------------------------

function phi_interp( d1, d2 )

! Finds soft phi(d1,d2) kernel by interpolation of phi_table

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: phi_interp

  integer :: i1, i2, id1, id2
  real(dp):: dd1(0:3), dd2(0:3)

  if (.not.phi_table_set) call set_phi_table()

  if (d1>=dcut .or. d2>=dcut) then
    phi_interp = 0
    return
  end if

  id1 = iofd( d1 )
  id2 = iofd( d2 )

  dd1(0) = 1
  dd2(0) = 1
  dd1(1) = (d1 - dmesh(id1)) / (dmesh(id1+1) - dmesh(id1))
  dd2(1) = (d2 - dmesh(id2)) / (dmesh(id2+1) - dmesh(id2))
  dd1(2) = dd1(1)**2
  dd2(2) = dd2(1)**2
  dd1(3) = dd1(1)**3
  dd2(3) = dd2(1)**3

  phi_interp = 0
  do i2 = 0,3
    do i1 = 0,3
      phi_interp = phi_interp + phi_table(i1,i2,id1,id2) * dd1(i1) * dd2(i2)
    end do
  end do

end function phi_interp

! -----------------------------------------------------------------------------

subroutine phiofr( r, phi )

! Finds phi(q1,q2,r) = phi(q1*r,q2*r) with q1=qmesh(iq1), q2=qmesh(iq2), 
! by interpolation from phi_table
! Notice that q is a density parameter, related to the Fermi wavevector

  implicit none
  real(dp),intent(in) :: r
  real(dp),intent(out):: phi(:,:)

  integer :: iq1, iq2
  real(dp):: dphidr

  if (size(phi,1)<mq .or. size(phi,2)<mq) &
    stop 'phiofr: ERROR: size(phi) too small'
  if (.not.qmesh_set) call set_qmesh()

  if (r >= rcut) then
    phi(:,:) = 0
  else
    do iq2 = 1,mq
      do iq1 = 1,iq2
        ! Use unfiltered kernel
!        phi(iq1,iq2) = phi_interp( qmesh(iq1)*r, qmesh(iq2)*r )
        ! Use filtered kernel
        call splint( dr, phir(:,iq1,iq2), d2phidr2(:,iq1,iq2), nr+1, r, &
                     phi(iq1,iq2), dphidr )
        phi(iq2,iq1) = phi(iq1,iq2)
      end do ! iq1
    end do ! iq2
  end if ! (r>=rcut)

end subroutine phiofr

!-----------------------------------------------------------------------------

function phi_soft( d1, d2 )

! Finds phi(d1,d2) softened near d1=d2=0

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: phi_soft

  real(dp):: d, d1m, d1p, d2m, d2p, dphidd, &
             phi, phi0, phi2, phi4, phim, phip

  d = sqrt( d1**2 + d2**2 )

  if (d<=0._dp) then
    phi_soft = phi0_soft
  else if (d > dsoft) then
    phi_soft = phi_val( d1, d2 )
  else ! (0<d<dsoft)

    d1p = d1/d * (dsoft+dd)
    d2p = d2/d * (dsoft+dd)
    d1m = d1/d * (dsoft-dd)
    d2m = d2/d * (dsoft-dd)
    phip = phi_val( d1p, d2p )
    phim = phi_val( d1m, d2m )
    phi = (phip + phim) / 2
    dphidd = (phip - phim) / (2*dd)
!    phi0 = phi - dphidd*dsoft/2
    phi0 = phi0_soft
    phi2 = (4*(phi-phi0) - dphidd*dsoft) / (2*dsoft**2)
    phi4 = (2*(phi0-phi) + dphidd*dsoft) / (2*dsoft**4)
    phi_soft = phi0 + phi2*d**2 + phi4*d**4
! BEGIN DEBUG
!    print'(a,5f8.3)', 'phi_soft: d,delta,phi0,phi2,phi4=', &
!      d, (d1-d2)/(d1+d2), phi0, phi2, phi4
! END DEBUG

  end if ! (d<=0)

end function phi_soft

! -----------------------------------------------------------------------------

real(dp) function phi_val( d1, d2 )

! Finds kernel phi by direct integration

  real(dp),intent(in) :: d1, d2

  integer  :: ia, ib, n, nmesh
  real(dp) :: a, acut, b, da1, dan, deff, dmax, dmin, gamma, &
              pi, t, w
  real(dp),allocatable :: amesh(:), c(:), dphida(:), dphidb(:), &
                          s(:), nu1(:), nu2(:)

! Find integration mesh
  dmax = max( abs(d1), abs(d2) )
  dmin = min( abs(d1), abs(d2) )
!  deff = dmax
  deff = sqrt(d1**2+d2**2)
  acut = max( acutmin, acutbyd * deff )
  da1 = min( damax, dabyd * deff )
  da1 = max( da1, damin )
  dan = damax
  n = get_n( 0._dp, acut, da1, dan )
  allocate( amesh(n), c(n), dphida(n), dphidb(n), s(n), nu1(n), nu2(n) )
  call set_mesh( n, xmax=acut, dxndx1=dan/da1 )
  call get_mesh( n, nmesh, amesh )

! BEGIN DEBUG
!  print'(a,i6,/,(10f8.3))', 'phi_val: size, amesh =', n, amesh
! END DEBUG

! Find cos(a), sin(a), nu1(a), and nu2(a)
  pi = acos(-1._dp)
  gamma = 4*pi/9
  do ia = 1,n
    a = amesh(ia)
    c(ia) = cos(a)
    s(ia) = sin(a)
    if (d1==0._dp) then
      nu1(ia) = a*a / 2
    else
      nu1(ia) = a*a / 2 / (1-exp(-gamma*(a/d1)**2))
    end if
    if (d2==0._dp) then
      nu2(ia) = a*a / 2
    else
      nu2(ia) = a*a / 2 / (1-exp(-gamma*(a/d2)**2))
    end if
  end do

! Make integral on variable a
  dphida(1) = 0
  do ia = 2,n
    a = amesh(ia)

    ! Make integral on variable b
    dphidb(1) = 0
    do ib = 2,n
      b = amesh(ib)

      w = 2*( (3-a*a)*b*c(ib)*s(ia) + (3-b*b)*a*c(ia)*s(ib) &
            + (a*a+b*b-3)*s(ia)*s(ib) - 3*a*b*c(ia)*c(ib) )/(a*b)**3

      t = 0.5_dp * ( 1/(nu1(ia)+nu1(ib)) + 1/(nu2(ia)+nu2(ib)) ) &
                 * ( 1/(nu1(ia)+nu2(ia))/(nu1(ib)+nu2(ib)) &
                   + 1/(nu1(ia)+nu2(ib))/(nu2(ia)+nu1(ib)) )

      dphidb(ib) = a*a * b*b * w * t

    end do ! ib
    dphida(ia) = integral( n, dphidb )

  end do ! ia
  phi_val = 2/pi**2 * integral( n, dphida )

  deallocate( amesh, c, dphida, dphidb, s, nu1, nu2 )

! BEGIN DEBUG
!  print'(a,2f12.6,i8,f12.6)', 'phi_val: d1,d2, na, phi =', d1, d2, n, phi_val
! END DEBUG

end function phi_val

! -----------------------------------------------------------------------------

subroutine pofq( q0, p0, dp0dq0 )

! Finds the values and derivatives, at q0, of the cubic polynomials 
! p_i(q0) such that
!    y(q0) = Sum_i p_i(q0) * y_i
! is the cubic spline interpolation at q0 of (any) function y(q) with
! values y_i at mesh points qmesh_i

  implicit none
  real(dp),intent(in) :: q0
  real(dp),intent(out):: p0(mq)
  real(dp),intent(out):: dp0dq0(mq)

  integer :: iq, iq0
  real(dp):: a, b, dq
  logical, save :: first_call=.true.
  real(dp),save :: p(mq,mq), d2pdq2(mq,mq)

! Set up spline polynomial basis
  if (first_call) then
    p = 0
    do iq = 1,mq
      p(iq,iq) = 1
      call generate_spline( qmesh, p(:,iq), mq, d2pdq2(:,iq) )
!      call generate_spline( qmesh, p(:,iq), mq, d2pdq2(:,iq), 0._dp, 0._dp )
    end do
    first_call = .false.
  end if

! Find interval of qmesh in which q0 is included
  if (q0>qmesh(mq)) then   ! q0 out of range
    p0 = 0
    dp0dq0 = 0
    return
  end if
  iq0 = iofq( q0 )

! Evaluate polynomials of spline basis
  dq = qmesh(iq0+1) - qmesh(iq0)
  a = (qmesh(iq0+1) - q0) / dq   ! dadq0 = -1/dq
  b = (q0 - qmesh(iq0)) / dq     ! dbdq0 = +1/dq
  do iq = 1,mq
    p0(iq) = a*p(iq0,iq) + b*p(iq0+1,iq) &
      + ((a**3-a)*d2pdq2(iq0,iq) + (b**3-b)*d2pdq2(iq0+1,iq)) * dq**2/6
    dp0dq0(iq) = - (p(iq0,iq) - p(iq0+1,iq)) / dq &
      - ((3*a**2-1)*d2pdq2(iq0,iq) - (3*b**2-1)*d2pdq2(iq0+1,iq)) * dq/6
  end do

end subroutine pofq

!-----------------------------------------------------------------------------

subroutine qofrho( rho, grho, q, dqdrho, dqdgrho )

! Finds the local wavevector parameter q0 defined in eqs.(11-12) of
! M.Dion et al, PRL 92, 246401 (2004)

  implicit none
  real(dp), intent(in) :: rho        ! Electron density
  real(dp), intent(in) :: grho(3)    ! Density gradient
  real(dp), intent(out):: q          ! Local wave vector parameter q0
  real(dp), intent(out):: dqdrho     ! d_q/d_rho
  real(dp), intent(out):: dqdgrho(3) ! d_q/d_grho

  character(len=*),parameter :: author = 'PW92'  ! Perdew-Wang'92 for LDA
  integer,         parameter :: irel   = 0       ! Non-relativistic exchange
  integer,         parameter :: nspin  = 1       ! Unpolarized electron gas
  real(dp),parameter :: zab = -0.8491_dp         ! See Dion et al
  real(dp):: decdrho, dexdrho, dkfdrho, dq0dgrho(3), dq0dgrho2, dq0drho, &
             dqdq0, dvxdrho, dvcdrho, ex, ec, grho2, kf, pi, q0, vx, vc

! Trap exception for zero density
  if (rho <= 1.e-15_dp) then
    q = qcut
    dqdrho = 0
    dqdgrho = 0
    return
  end if

  pi = acos(-1._dp)
  kf = (3*pi**2 * rho)**(1._dp/3)

! Find exchange and correlation energy densities
  call ldaxc( author, irel, nspin, rho, ex, ec, vx, vc, dvxdrho, dvcdrho )

! Find q
  grho2 = sum(grho**2)
  q0 = ( 1 + ec/ex - zab/9 * grho2 / (2*kf*rho)**2 ) * kf

! Find derivatives
  dkfdrho = kf / (3*rho)
  dexdrho = (vx - ex) / rho  ! Since vx = d(rho*ex)/d_rho = ex + rho*dexdrho
  decdrho = (vc - ec) / rho
  dq0drho = ( decdrho/ex - ec/ex**2*dexdrho + 2 * zab/9 * grho2 / &
              (2*kf*rho)**3 * (2*dkfdrho*rho + 2*kf) ) * kf &
            + q0/kf * dkfdrho
  dq0dgrho2 = -(zab/9) / (2*kf*rho)**2 * kf
  dq0dgrho(:) = dq0dgrho2 * 2*grho(:)  ! Since d(vector**2)/d_vector = 2*vector

! Saturate q to qcut smoothly
  call saturate( q0, qcut, q, dqdq0 )
  dqdrho = dqdq0 * dq0drho
  dqdgrho = dqdq0 * dq0dgrho

end subroutine qofrho

!-----------------------------------------------------------------------------

subroutine saturate( x, xc, y, dydx )

  ! Defines a function y(x) = xc * (1 - exp(-Sum_n=1:nsat (x/xc)**n/n))
  ! It is approx. equal to x for x<xc and it saturates to xc when x->infinity

  implicit none
  real(dp),intent(in) :: x     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: y     ! Function value
  real(dp),intent(out):: dydx  ! Derivative dy/dx

  integer :: n
  real(dp):: dpdx, p

  if (nsat >= 100) then
    if (x < xc) then
      y = x
      dydx = 1
    else ! (x >= xc)
      y = xc
      dydx = 0
    end if ! (x < xc)
  else ! (nsat < 100)
!    This is the straightforward polynomial evaluation
!    p = x/xc
!    dpdx = 1/xc
!    do n = 2,nsat
!      p = p + (x/xc)**n / n
!      dpdx = dpdx + (x/xc)**(n-1) / xc
!    end do
!   And this is a more accurate way to evaluate it
    p = (x/xc)/nsat
    dpdx = 1._dp/xc
    do n = nsat-1,1,-1
      p = (p + 1._dp/n) * x/xc
      dpdx = (dpdx*x + 1) / xc
    end do
    y = xc * (1 - exp(-p))
    dydx = xc * dpdx * exp(-p)
  end if ! (nsat >= 100)

end subroutine saturate

!-----------------------------------------------------------------------------

subroutine saturate_inverse( y, xc, x, dydx )

! Finds the inverse of saturate function

  implicit none
  real(dp),intent(in) :: y     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: x     ! Inverse function value
  real(dp),intent(out):: dydx  ! Derivative dy/dx

  real(dp):: x1, x2, yx

  if (y<0._dp .or. y>xc) stop 'vdw:saturate_inverse: y out of range'
  x1 = 0
  x2 = xmaxbyxc * xc
  do
    x = (x1+x2)/2
    call saturate( x, xc, yx, dydx )
    if (abs(y-yx)<ytol) then
      return
    else if (yx < y) then
      x1 = x
    else
      x2 = x
    end if
  end do

end subroutine saturate_inverse

!-----------------------------------------------------------------------------

subroutine set_phi_table()

! Finds the interpolation table (mesh points and function values) for phi(d1,d2)

  implicit none

  logical :: file_found
  integer :: id, id1, id2, nmesh
  real(dp):: d, d1, d1m, d1p, d2, d2m, d2p, dd, &
             dphidd1(nd,nd), dphidd2(nd,nd), d2phidd1dd2(nd,nd), &
             phi(nd,nd), phi1, phim, phip, phimm, phimp, phipm, phipp

! Read file with table, if present
  inquire( file='phi.table', exist=file_found )
  if (file_found) then
    open( unit=1, file='phi.table', form='unformatted' )
    read(1,end=1) nmesh
    if (nmesh==nd) then
      read(1,end=1) dmesh
      read(1,end=1) phi_table
    end if
    close( unit=1 )
    phi_table_set = .true.
    return
  end if
1 continue ! Come here if end-of-file found

! Set d-mesh points
  call set_mesh( nd, xmax=dcut, dxndx1=ddmaxddmin )
  call get_mesh( nd, nmesh, dmesh )
! BEGIN DEBUG
  write(udebug,'(a,/,(10f8.4))') 'set_phi_table: dmesh =', dmesh
! BEGIN DEBUG

! Find function at mesh points
  do id1 = 1,nd
    d1 = dmesh(id1)
    phi1 = phi_soft( d1, d1 )
    phi(id1,id1) = phi1
    do id2 = 1,id1-1
      d2 = dmesh(id2)
      d = sqrt( d1**2 + d2**2 )
      if (d < dsoft) then
        phi(id1,id2) = phi_soft( d1, d2 )
      else
        if (use_dphi) then ! Use dphi for better efficiency
          phi(id1,id2) = phi1 + dphi( d1, d2 )
        else ! Use only phi_val, to eliminate uncertainties
          phi(id1,id2) = phi_val( d1, d2 )
        end if
      end if
      phi(id2,id1) = phi(id1,id2)
    end do ! id2
  end do ! id1

  open( unit=4, file='phi.out' )
  do id2 = 1,nd
    write(4,'(/,(f12.6))') phi(:,id2)
  end do
  close( unit=4 )

  if (deriv_method == 'numeric') then
!    print*, 'set_phi_table: Using numerical derivatives'

    ! Find derivatives at mesh points
     do id1 = 1,nd
      d1 = dmesh(id1)
      dd = ddbyd * d1
      dd = max( dd, ddmin )
!      d1 = max( d1, dd )
      d1m = d1 - dd
      d1p = d1 + dd
      phim = phi_soft( d1m, d1m )
      phip = phi_soft( d1p, d1p )
      do id2 = 1,id1
        d2  = dmesh(id2)
!        d2 = max( d2, dd )
        d = sqrt( d1**2 + d2**2 )
        d1m = d1 - dd
        d1p = d1 + dd
        d2m = d2 - dd
        d2p = d2 + dd

        if (d < dsoft) then
          phimm = phi_soft( d1m, d2m )
          phipm = phi_soft( d1p, d2m )
          phipp = phi_soft( d1p, d2p )
          phimp = phi_soft( d1m, d2p )
        else ! (d>dsoft)
          if (use_dphi) then
            phimm = phim + dphi( d1m, d2m )
            phipm = phip + dphi( d1p, d2m )
            phipp = phip + dphi( d1p, d2p )
            if (id1==id2) then
              phimp = phip + dphi( d1m, d2p )
            else
              phimp = phim + dphi( d1m, d2p )
            end if
          else ! (.not.use_dphi)
            phimm = phi_val( d1m, d2m )
            phipm = phi_val( d1p, d2m )
            phipp = phi_val( d1p, d2p )
            phimp = phi_val( d1m, d2p )
          end if ! (use_dphi)
        end if ! (d<dsoft)

        dphidd1(id1,id2)     = (phipp+phipm-phimp-phimm) / (4*dd)
        dphidd2(id1,id2)     = (phipp-phipm+phimp-phimm) / (4*dd)
        d2phidd1dd2(id1,id2) = (phipp-phipm-phimp+phimm) / (2*dd)**2
        dphidd1(id2,id1)     = dphidd2(id1,id2)
        dphidd2(id2,id1)     = dphidd1(id1,id2)
        d2phidd1dd2(id2,id1) = d2phidd1dd2(id1,id2)

      end do ! id2
    end do ! id1

  else if (deriv_method == 'interp') then
!    print*, 'set_phi_table: Using interpolation for derivatives'

    ! Reset mesh, which has been changed by phi_val
    call set_mesh( nd, xmax=dcut, dxndx1=ddmaxddmin )

    ! Set interpolation method
    if (interp_method=='Lagrange') then
      call set_interpolation( 'Lagrange' )
    else if (interp_method=='Spline') then
!      call set_interpolation( 'Spline', huge(phi), huge(phi) )
      call set_interpolation( 'Spline', 0._dp, 0._dp )
    else
      stop 'set_phi_val: ERROR: Unknown interp_method'
    end if

    ! Find first partial derivative d_phi/d_d1
    do id = 1,nd
      dphidd1(:,id) = derivative( nd, phi(:,id) )
      dphidd2(id,:) = dphidd1(:,id)
    end do

    ! Find second cross partial derivative d_phi/d_d1/d_d2
    do id = 1,nd
      d2phidd1dd2(id,:) = derivative( nd, dphidd1(id,:) )
    end do

    ! Symmetrize d_phi/d_d1/d_d2
    do id2 = 2,nd
      do id1 = 1,id2-1
        d2phidd1dd2(id1,id2) = (d2phidd1dd2(id1,id2) + &
                                d2phidd1dd2(id2,id1)) / 2
        d2phidd1dd2(id2,id1) = d2phidd1dd2(id1,id2) 
      end do
    end do

  else
    stop 'set_phi_table: ERROR: Unknown deriv_method'
  end if ! (deriv_method)

! Make values and derivatives strictly zero when d1=dmax or d2=dmax
  phi(:,nd) = 0
  phi(nd,:) = 0
  dphidd1(:,nd) = 0
  dphidd1(nd,:) = 0
  dphidd2(:,nd) = 0
  dphidd2(nd,:) = 0
  d2phidd1dd2(:,nd) = 0
  d2phidd1dd2(nd,:) = 0

! Make dphi(d1,d2)/dd1=0 for d1=0 and dphi(d1,d2)/dd2=0 for d2=0
  dphidd1(1,:) = 0
  dphidd2(:,1) = 0
  d2phidd1dd2(:,1) = 0
  d2phidd1dd2(1,:) = 0

! BEGIN DEBUG
! Print values and derivatives for debugging
!  print'(a,/,2a10,4a15)', &
!   'set_phi_table:', 'd1', 'd2', 'phi', 'dphi/dd1', 'dphi/dd2', 'd2phi/dd1dd2'
!  do id1 = 1,nd
!    do id2 = id1,nd
!      d1 = dmesh(id1)
!      d2 = dmesh(id2)
!      print'(2f10.6,4e15.6)', d1, d2, phi(id1,id2), &
!        dphidd1(id1,id2), dphidd2(id1,id2), d2phidd1dd2(id1,id2)
!    end do
!  end do
! END DEBUG

! Set up bicubic interpolation coefficients
  call bcucof( nd, nd, dmesh, dmesh, phi, dphidd1, dphidd2, d2phidd1dd2, &
               phi_table )

! Save phi_table in file
  open( unit=1, file='phi.table', form='unformatted' )
  write(1) nd
  write(1) dmesh
  write(1) phi_table
  close( unit=1 )

! Mark table as set
  phi_table_set = .true.

end subroutine set_phi_table

! -----------------------------------------------------------------------------

subroutine set_qmesh()

! Sets mesh of q values

  implicit none
  integer :: nmesh

! BEGIN DEBUG
!      real(dp):: a, b
!      a = log(20.0_dp) / (20-1)
!      b = 8.0_dp / (exp(a*(20-1)) - 1)
!      qcut = b * (exp(a*(nq-1)) - 1)
!      dqmaxdqmin = exp(a*(nq-1))
! END DEBUG

  call set_mesh( mq, xmax=qcut, dxndx1=dqmaxdqmin )
  call get_mesh( mq, nmesh, qmesh )
  qmesh_set = .true.

! BEGIN DEBUG
  write(udebug,'(/,a,/,(10f8.4))') 'vdw:set_qmesh: qmesh =', qmesh
! END DEBUG

end subroutine set_qmesh

! -----------------------------------------------------------------------------

subroutine vdw_decusp( nspin, rhos, grhos, eps, dedrho, dedgrho )

  implicit none
  integer, intent(in) :: nspin            ! Number of spin components
  real(dp),intent(in) :: rhos(nspin)      ! Electron spin density
  real(dp),intent(in) :: grhos(3,nspin)   ! Spin density gradient
  real(dp),intent(out):: eps              ! Energy correction per electron
  real(dp),intent(out):: dedrho(nspin)    ! d(rho*eps)/d(rho)
  real(dp),intent(out):: dedgrho(3,nspin) ! d(rho*eps)/d(grad_rho)

  logical, save:: initialized=.false.
  real(dp),save:: table(mq,mq)
  integer :: iq1, iq2, ir, is, ix, ns
  real(dp):: dptpdq, dpdq(mq), dqdrho, dqdgrho(3), grho(3), p(mq), &
             ptp, phi22, phi(mq,mq), pi, pt(mq), q, r, rho

  if (.not.initialized) then

    if (.not.phi_table_set) call set_phi_table()
    if (.not.kcut_set) stop 'vdw_decusp: ERROR: kcut has not been set'

    pi = acos(-1._dp)
    table = 0
    do ir = 1,nr
      r = ir * dr
      call phiofr( r, phi )
      do iq2 = 1,mq
        do iq1 = 1,iq2
          table(iq1,iq2) = table(iq1,iq2) - 2*pi*dr * r**2 * phi(iq1,iq2)
        end do ! iq1
      end do ! iq2
    end do

    do iq2 = 2,mq
      do iq1 = 1,iq2-1
        table(iq2,iq1) = table(iq1,iq2)
      end do
    end do

    initialized = .true.

  end if ! (.not.initialized)

  ns = min(nspin,2)
  rho = sum(rhos(1:ns))
  do ix = 1,3
    grho(ix) = sum(grhos(ix,1:ns))
  end do

  call qofrho( rho, grho, q, dqdrho, dqdgrho )
  call pofq( q, p, dpdq )

  pt = matmul(p,table)
  ptp = sum( pt * p )
  dptpdq = 2 * sum( pt * dpdq )
  eps = rho * ptp
  dedrho(:) = 2 * rho * ptp + rho**2 * dptpdq * dqdrho
  do ix = 1,3
    dedgrho(ix,:) = rho**2 * dptpdq * dqdgrho(ix)
  end do

end subroutine vdw_decusp

!-----------------------------------------------------------------------------

subroutine vdw_get_qmesh( n, q )

! Returns size and values of q-mesh

  implicit none
  integer,          intent(out) :: n
  real(dp),optional,intent(out) :: q(:)
  integer:: nmax
  if (.not.qmesh_set) call set_qmesh()
  n = nq
  if (present(q)) then
    nmax = max( nq, size(q) )
    q(1:nmax) = qmesh(1:nmax)
  end if
end subroutine vdw_get_qmesh

! -----------------------------------------------------------------------------

subroutine vdw_phi( k, phi )

! Finds and interpolates phi(q1,q2,k) (Fourier transform of phi(q1,q2,r)) 
! for all values of q1 and q2 in qmesh, in the first call. Then it finds
! the interpolation at that and succesive calls

  implicit none
  real(dp),intent(in) :: k         ! Modulus of actual k vector
  real(dp),intent(out):: phi(:,:)  ! phi(q1,q2,k) at given k
                                   ! for all q1,q2 in qmesh

  integer :: ik, iq1, iq2
  real(dp):: a, a3, b, b3, dphidk

  if (.not.kcut_set) stop 'vdw_phi: ERROR: kcut must be previously set'

! Check argument sizes
  if (size(phi,1)<nq .or. size(phi,2)<nq) &
    stop 'vdw_phi: ERROR: size(phi) too small'

! Find phi values at point k
  if (k >= kcut) then
    phi(:,:) = 0
  else
    ! Expand interpolation inline since this is the hottest point in VdW
    ik = k/dk
    a = ((ik+1)*dk-k)/dk
    b = 1 - a
    a3 = (a**3 - a) * dk**2 / 6
    b3 = (b**3 - b) * dk**2 / 6
    do iq2 = 1,nq
      do iq1 = 1,iq2
!        call splint( dk, phik(:,iq1,iq2), d2phidk2(:,iq1,iq2), nk+1, k, &
!                     phi(iq1,iq2), dphidk, pr )
        phi(iq1,iq2) = a*phik(ik,iq1,iq2) + b*phik(ik+1,iq1,iq2) &
                + a3*d2phidk2(ik,iq1,iq2) + b3*d2phidk2(ik+1,iq1,iq2)
        phi(iq2,iq1) = phi(iq1,iq2)
      end do
    end do
  end if

end subroutine vdw_phi

!-----------------------------------------------------------------------------

subroutine vdw_set_kcut( kc )

! Sets the reciprocal planewave cutoff kc of the integration grid, and finds 
! the interpolation table to be used by vdw_phi to obtain the vdW kernel phi
! at the reciprocal mesh wavevectors.

  implicit none
  real(dp),intent(in):: kc  ! Planewave cutoff: k>kcut => phi=0

  integer :: ik, iq1, iq2, ir, nrs
  real(dp):: dphids, dphidk0, dphidkmax, dphidr0, dphidrmax, dqdq0, &
             k, kmax, phi(0:nr), phi0, phi2, phis, pi, q1, q2, r(0:nr), rs

  if (kc == kcut) return   ! Work alredy done
  if (.not.qmesh_set) call set_qmesh()

  pi = acos(-1._dp)
  dr = rcut / nr
  dk = pi / rcut
  kmax = pi / dr
  nk = int(kc/dk) + 1
  nrs = nint(rsoft/dr)
  rs = nrs * dr

! BEGIN DEBUG
  write(udebug,'(a,5f8.3)') 'vdw_set_kcut: dcut,qcut,rcut,kcut,kmax=', &
    dcut, qcut, rcut, kc, kmax
! END DEBUG

  ! For each pair of values q1 and q2
  do iq2 = 1,mq
    do iq1 = 1,iq2
!      print*, 'vdw_set_kcut: iq1,iq2=', iq1, iq2

!     Saturated q values
!      q1 = qmesh(iq1)
!      q2 = qmesh(iq2)

!     Find original (unsaturated) q values
      call saturate_inverse( qmesh(iq1), qcut, q1, dqdq0 )
      call saturate_inverse( qmesh(iq2), qcut, q2, dqdq0 )

      ! Find kernel in real space
      do ir = 0,nr
        r(ir) = ir * dr
        phir(ir,iq1,iq2) = phi_interp( q1*r(ir), q2*r(ir) )
      end do
      phi(:) = phir(:,iq1,iq2)

      ! Change kernel near origin to a parabola matching at rs
      if (nrs>0) then
        phis = phi(nrs)
        dphids = (phi(nrs+1) - phi(nrs-1)) / (2*dr)
        phi0 = phis - dphids * rs/2
        phi2 = dphids / (2*rs)
        do ir = 0,nrs
          phir(ir,iq1,iq2) = phi0 + phi2 * r(ir)**2
        end do
      end if ! (nrs>0)

      ! Kill kernel smoothly at r=rcut
      do ir = 0,nr
        phir(ir,iq1,iq2) = phir(ir,iq1,iq2) * cutoff( r(ir)/rcut )
      end do

      ! Optimized filter (warning: inaccurate for very large kcut*rcut)
!      call filter( 0, nr+1, r(:), phir(:,iq1,iq2), kc, 0 )

      ! Find kernel in reciprocal space
      call radfft( 0, nr, rcut, phir(:,iq1,iq2), phik(:,iq1,iq2) )
      phik(:,iq1,iq2) = phik(:,iq1,iq2) * (2*pi)**1.5_dp

      ! Filter out above kcut
      phik(nk:nr,iq1,iq2) = 0

      ! Soft filter below kcut
      do ik = 1,nk
        k = ik * dk
        phik(ik,iq1,iq2) = phik(ik,iq1,iq2) * cutoff(k/kc)
      end do

      ! Find filtered kernel in real space
      call radfft( 0, nr, kmax, phik(:,iq1,iq2), phir(:,iq1,iq2) )
      phir(:,iq1,iq2) = phir(:,iq1,iq2) / (2*pi)**1.5_dp

      ! Set up spline interpolation tables
      dphidr0 = 0
      dphidrmax = 0
      dphidk0 = 0
      dphidkmax = 0
      call spline( dr, phir(:,iq1,iq2), nr+1, dphidr0, dphidrmax, &
                   d2phidr2(:,iq1,iq2) )
      call spline( dk, phik(:,iq1,iq2), nk+1, dphidk0, dphidkmax, &
                   d2phidk2(:,iq1,iq2) )

      ! Fill symmetric elements
      phir(:,iq2,iq1) = phir(:,iq1,iq2)
      phik(:,iq2,iq1) = phik(:,iq1,iq2)
      d2phidr2(:,iq2,iq1) = d2phidr2(:,iq1,iq2)
      d2phidk2(:,iq2,iq1) = d2phidk2(:,iq1,iq2)

!      if (.false. .and. iq1==iq2) then
!        print*, 'vdw_set_kcut: iq1,iq2=', iq1, iq2
!        call window( 0._dp, 5._dp, -1._dp, 4._dp, 0 )
!        call axes( 0._dp, 1._dp, 0._dp, 1._dp )
!        call plot( nr+1, r, phi, phir(:,iq1,iq2) )
!        call window( 0._dp, 10._dp, -0.05_dp, 0.15_dp, 0 )
!        call axes( 0._dp, 1._dp, 0._dp, 0.05_dp )
!        call plot( nr+1, r, q1*q2*r**2*phi, q1*q2*r**2*phir(:,iq1,iq2) )
!        call show()
!      end if

    end do ! iq1
  end do ! iq2

!  print'(a,/,(2i6,f12.6))', 'vdw_set_kcut: iq1, iq2, phir(0,iq1,iq2) =', &
!    ((iq1,iq2,phir(0,iq1,iq2),iq1=2,iq2),iq2=2,mq)

  kcut = kc
  kcut_set = .true.

end subroutine vdw_set_kcut

!-----------------------------------------------------------------------------

subroutine vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )

  implicit none
  integer, intent(in) :: nspin               ! Number of spin components
  real(dp),intent(in) :: rhos(nspin)         ! Electron spin density
  real(dp),intent(in) :: grhos(3,nspin)      ! Spin density gradient
  real(dp),intent(out):: theta(nq)           ! Expansion of rho*q in qmesh
  real(dp),intent(out):: dtdrho(nq,nspin)    ! dtheta(iq)/drhos
  real(dp),intent(out):: dtdgrho(3,nq,nspin) ! dtheta(iq)/dgrhos(ix)

  integer :: is, ix, ns
  real(dp):: rho, grho(3), dpdq(mq), dqdrho, dqdgrho(3), p(mq), q

  ns = min(nspin,2)
  rho = sum(rhos(1:ns))
  do ix = 1,3
    grho(ix) = sum(grhos(ix,1:ns))
  end do

  call qofrho( rho, grho, q, dqdrho, dqdgrho )
  call pofq( q, p, dpdq )

  theta(1:nq) = rho * p(1:nq)
  dtdrho(:,:) = 0
  dtdgrho(:,:,:) = 0
  do is = 1,ns
    dtdrho(1:nq,is) = p(1:nq) + rho * dpdq(1:nq) * dqdrho
    do ix = 1,3
      dtdgrho(ix,1:nq,is) = rho * dpdq(1:nq) * dqdgrho(ix)
    end do
  end do

end subroutine vdw_theta

END MODULE m_vdwxc
