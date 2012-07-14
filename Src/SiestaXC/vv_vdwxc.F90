!!@LICENSE
!
!******************************************************************************
! MODULE m_vv_vdwxc
! Implements the nonlocal correlation energy part of the van der Waals density
! functional of O.A.Vydrov & T. van Voorhis, JCP 133, 244103 (2010):
!   Enlc = (1/2) Int Int dr1 dr2 n(r1) phi(n1,gn1,n2,gn2,r12) n(r2)
! where r12=|r2-r1|, n(r) is the electron density at point r, gn(r) its 
! gradient, and phi is a universal function defined in Eqs.(2-6) and (9).
! To be used with module m_vdwxc
! Refs: O.A.Vydrov & T.vanVoorhis, JCP 133, 244103 (2010)
!       G.Roman-Perez and J.M.Soler, PRL 103, 096102 (2009)
! Written by J.M.Soler. July 2012
!------------------------------------------------------------------------------
! Used module procedures:
!  use flib_spline, only: generate_spline   ! Sets spline in a general mesh
!  use mesh1D,      only: get_mesh          ! Returns the mesh points
!  use m_radfft,    only: radfft            ! Radial fast Fourier transform
!  use mesh1D,      only: set_mesh          ! Sets a 1D mesh
!  use m_recipes,   only: spline            ! Sets spline in a uniform mesh
!  use m_recipes,   only: splint            ! Performs spline interpolation
!------------------------------------------------------------------------------
! Used module parameters:
!   use precision,  only: dp                ! Real double precision type
!------------------------------------------------------------------------------
! Public procedures available from this module:
!   vv_vdw_beta      : Returns constant beta of eq.(13) of VV2010 JCP paper
!   vv_vdw_get_kmesh : Returns size and values of (kf,kg) mesh
!   vv_vdw_phi       : Finds phi(k) at (kf,kg) mesh points
!   vv_vdw_set_kcut  : Sets the planewave cutoff kc of the integration grid
!   vv_vdw_theta     : Finds function theta(n,gn) (eq.(8) of Roman-Soler PRL)
!------------------------------------------------------------------------------
! Public types, parameters, variables, and arrays:
!   None
!------------------------------------------------------------------------------
! Units: 
!   All lengths and energies in atomic units (Bohr, Hartree)
!******************************************************************************
! real(dp) function vv_vdw_beta()
!   Returns parameter beta=0.00497 of eq.(13) of VV2010 JCP paper
! Arguments:
!   none
! Sample usage:
!   real(dp):: beta
!   beta = vv_vdw_beta()
!------------------------------------------------------------------------------
! subroutine vdw_get_kmesh( nkf, nkg, kf, kg )
!   Returns size and values of the (kf,kg) interpolation mesh of the
!   Roman-Soler method. kf and kg are local Fermi- and gradient-wavevectors:
!   kf=(3*pi**2*n)**(1/3), kg=|grad(n)|/n
! Arguments:
!   integer,          intent(out) :: nkf    ! Number of kf mesh points
!   integer,          intent(out) :: nkg    ! Number of kg mesh points
!   real(dp),optional,intent(out) :: kf(:)  ! Mesh values of Fermi wavevect.
!   real(dp),optional,intent(out) :: kg(:)  ! Mesh values of grad(n)/n
! Sample usage:
!   integer :: nkf, nkg
!   real(dp):: kcut
!   real(dp),allocatable:: kfmesh(:), kgmesh(:)
!   call vv_vdw_get_kmesh( nkf, nkg )
!   allocate( kfmesh(nkf), kgmesh(nkg) )
!   call vv_vdw_get_kmesh( nkf, nkg, kfmesh, kgmesh )
! Notes:
! - If the size arrays kf or kg is smaller than that of their stored meshes, 
!   they are filled with the first size(kf) or size(kg) values
! - The size and values of the logarithmic k meshes are set by internal
!   parameters that can be changed only by editing them in this module:
!     nkf, nkg          : Number of kf and kg mesh points
!     kfcut=kfmesh(nkf) : Max. value of kf mesh
!     kgcut=kgmesh(nkg) : Max. value of kg mesh
!     dkmaxdkmin     : (kmesh(nk) - kmesh(nk-1)) / (kmesh(2) - kmesh(1))
!                      (applicable to both kf and kg meshes)
!   Although the presently-set values have been found to yield good accuracy
!   in preliminary calculations, more tests may be required to guarantee
!   convergence in other systems. Larger nk's increase accuracy but CPU time 
!   increases between (nkf*nkg) and (nkf*nkg)**2.
!------------------------------------------------------------------------------
! subroutine vv_vdw_phi( k, phi, dphidk )
!   Finds phi(k1,k2,k) (Fourier transform of phi_soft(k1,k2,r)) for all 
!   values of k1=(kf1,kg1) and k2=(kf2,kg2) of the (kf,kg) mesh. The mesh
!   indexes of kf and kg are combined into a single index ikfg=1,...,nkf*nkg
!   so that the size of phi and dphidk is (nkf*nkg,nkg*nkg).
!   In practice, it returns n1*n2*phi(k1,k2,k), where n1=kf1**3/(3*pi**2) and
!   n1=kf1**3/(3*pi**2) at the mesh values of kf1 and kf2. This is because
!   this function is smoother to interpolate than phi itself.
! Arguments:
!   real(dp),intent(in) :: k            ! Modulus of actual k vector
!   real(dp),intent(out):: phi(:,:)     ! phi(k1,k2,k) at given k
!                                       ! for all k1,k2 in kmesh
!   real(dp),intent(out):: dphidk(:,:)  ! dphi(k1,k2,k)/dk at given k
! Sample usage:
!   integer :: nkf, nkg
!   real(dp):: k, kcut
!   real(dp),allocatable:: phi(:,:), dphidk(:,:)
!   kcut = 10._dp ! 10 Bohr^-1 => 100 Ry (this should be adapted to your mesh)
!   call vv_vdw_set_kcut( kcut )
!   call vv_vdw_get_qmesh( nkf, nkg )
!   allocate( phi(nkf*nkg,nkf*nkg), dphidk(nkf*nkg,nkf*nkg) )
!   do k points
!     call vv_vdw_phi( k, phi, dphidk )
! Notes:
! - Requires a previous call to vv_vdw_set_kcut. Otherwise stops with an error.
! - Stops with an error message if size of array phi is smaller than
!   (nkf*nkg)**2.
!-----------------------------------------------------------------------------
! subroutine vv_vdw_set_kcut( kc )
!   Sets the reciprocal planewave cutoff kc of the integration grid, and finds 
!   the interpolation table to be used by vv_vdw_phi to obtain the vdW kernel 
!   phi at the reciprocal mesh wavevectors.
! Arguments:
!   real(dp),intent(in):: kc  ! Planewave cutoff: k>kcut => phi=0
!------------------------------------------------------------------------------
! subroutine vv_vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )
!   Finds the value and derivatives of function theta_ikfg(rho,grad_rho) 
!   of eq.(8) of the Roman-Soler PRL. ikfg is a combined index for the (kf,kg) 
!   interpolation mesh. 
!   In practice, it returns p_ikfg, rather than theta_ikgf=rho*p_ikgf, while
!   vv_vdw_phi returns rho1*rho2*phi, which is smoother to interpolate.
! Arguments:
!   integer, intent(in) :: nspin               ! Number of spin components
!   real(dp),intent(in) :: rhos(nspin)         ! Electron spin density
!   real(dp),intent(in) :: grhos(3,nspin)      ! Spin density gradient
!   real(dp),intent(out):: theta(nkf*nkg)      ! Function theta at mesh points
!   real(dp),intent(out):: dtdrho(nkf*nkg,nspin)     ! dtheta/drhos
!   real(dp),intent(out):: dtdgrho(3,nnkf*nkg,nspin) ! dtheta/dgrhos
! Notes:
! - Requires a previous call to vv_vdw_set_kcut
! - The values of kf(rho,grad_rho) and kg(rho,grad_rho) are saturated at
!   kfcut and kgcut (two parameters) to avoid incorrect 'interpolation' of 
!   phi(k1,k2,r12) from (kf,kg) mesh points.
!******************************************************************************
! Algorithms:
!   Once we set the cutoff kcut (expressed as a wavevector), of the integration 
! mesh to be used in evaluating the VdW nonlocal energy Enlc, we calculate and
! store (in memory) an interpolation table of phi_soft(k1,k2,k) for all values 
! of k1=(kf1,kg1) and k2=(kf2,kg2) in the (kf,kg) mesh, and of k in a fine 
! radial grid of nk=nr points. This is done by first calculating phi(k1,k2,r) 
! in a radial grid in real space and then Fourier-transforming it to k space.
! This table is then used to interpolate phi(k1,k2,k) for any desired k.
! In practice, we interpolate (and return) rho1*rho2*phi(k1,k2,k) (which is 
! smoother than phi(k1,k2,k) itself), where rho=kf**3/(3*pi**2).
!   In order to ensure that values (kf,kg) are within the interpolation range,
! they are 'saturated' smoothly to a cutoffs kfcut and kgcut. This implies an 
! approximation when either kf is very large (i.e. near the nucleus) or when 
! kg is very large, what tipically occurs in the tails of the electron density. 
! In the first case, Enlc is neglegible compared with Ex and the local part of 
! Ec. In the second case, it is neglegible because of the factor rho in the 
! integrand. Thus, the preliminary tests suggest that the presently-set values 
! of kfcut and kgcut give sufficiently accurate values.
!******************************************************************************

MODULE m_vv_vdwxc

! Used module procedures
  use sys,         only: die               ! termination routine
  use flib_spline, only: generate_spline   ! Sets spline in a general mesh
  use mesh1D,      only: get_mesh          ! Returns the mesh points
  use m_radfft,    only: radfft            ! Radial fast Fourier transform
  use mesh1D,      only: set_mesh          ! Sets a 1D mesh
  use m_recipes,   only: spline            ! Sets spline in a uniform mesh
  use m_recipes,   only: splint            ! Performs spline interpolation

! Used module parameters
  use precision,   only: dp                ! Real double precision type

#ifdef DEBUG_XC
  use debugXC,     only: udebug            ! File unit for debug output
!  use plot_module, only: plot
#endif /* DEBUG_XC */

  implicit none

! Called by m_vdwxc
PUBLIC ::            &
  vv_vdw_beta,       &! Returns parameter beta of the VV2010 functional
  vv_vdw_theta,      &! Finds function theta_ikfg(rho,grad_rho)
  vv_vdw_get_kmesh,  &! Returns size and values of (kf,kg) mesh
  vv_vdw_phi,        &! Finds and interpolates rho1*rho2*phi(k1,k2,k)
  vv_vdw_set_kcut     ! Sets the planewave cutoff kc of the integration grid

#ifdef DEBUG_XC
! Called by debugging test programs
PUBLIC :: &
  vv_vdw_phiofr       ! Kernel phi(k1,k2,r) at tabulated k1,k2-mesh values
#endif /* DEBUG_XC */

PRIVATE  ! Nothing is declared public beyond this point

!  integer, parameter:: dp = kind(1.d0)

  ! Set parameters of the Vydrov-vanVoorhis functional,
  ! given in paragraph after eq.(27) of JCP 2010 paper
  real(dp),parameter:: vv_C = 0.0093
  real(dp),parameter:: vv_b = 5.9
  real(dp),parameter:: vv_beta = 0.00497

  ! Set derivation methods to use for interpolation table
  character(len=*),parameter:: deriv_method = 'numeric'  !('numeric'|'interp')
  character(len=*),parameter:: interp_method= 'Spline' !('Lagrange'|'Spline')

  ! Mesh parameters for table of phi(k1,k2,r) and its Fourier transform
  integer, parameter:: nr = 2048             ! Radial mesh points (power of 2)
  real(dp),parameter:: rcut = 100._dp        ! Radial cutoff: r>rcut => phi=0
  real(dp),parameter:: rmin = 1.e-6_dp       ! Min. radius as denominator
  integer, parameter:: nkf = 6               ! Number of Fermi wavevectors
  integer, parameter:: nkg = 5               ! Num. of grad(n)/n wavevectors
  integer, parameter:: nkfg = nkf*nkg        ! Num. of (kf,kg) mesh points
  real(dp),parameter:: kfcut = 5.0_dp        ! Max. Fermi wavevec.
  real(dp),parameter:: kgcut = 5.0_dp        ! Max. grad(n)/n
  real(dp),parameter:: dkmaxdkmin = 20.0_dp  ! Last/first k mesh interval

  ! Parameters for cutoff function, used in radial Fourier transforms of phi
  integer, parameter:: ncut1 =  8      ! cutoff(x)=(1-x**ncut1)**ncut2
  integer, parameter:: ncut2 =  4

  ! Parameters for saturate function, used to enforce that q<qcut
  integer, parameter:: nsat = 12 ! xsat(x,xc)=1/(1/x**nsat+1/xc**nsat)**(1/nsat)

  ! Private module variables and arrays
  real(dp),save:: kfmesh(nkf)              ! Mesh points for Fermi wavevector
  real(dp),save:: kgmesh(nkg)              ! Mesh points for grad(n)/n
  logical, save:: phi_table_set=.false.    ! Has phi_table been set?
  logical, save:: kmesh_set=.false.        ! Has (kf,kg) mesh been set?
  logical, save:: kcut_set=.false.         ! Has kcut been set?
  real(dp),save:: phir(0:nr,nkfg,nkfg)     ! Table of phi(k1,k2,r)
  real(dp),save:: d2phidr2(0:nr,nkfg,nkfg) ! Table of d2_phi/dr2
  real(dp),save:: dr                       ! r-mesh interval
  real(dp),save:: phik(0:nr,nkfg,nkfg)     ! Table of phi(k1,k2,k)
  real(dp),save:: d2phidk2(0:nr,nkfg,nkfg) ! Table of d2_phi/dk2
  real(dp),save:: dk                       ! k-mesh interval
  real(dp),save:: kcut                     ! Planewave cutoff: k>kcut => phi=0
  real(dp),save:: kmax                     ! Max. k vector in Fourier transforms
  integer, save:: nk                       ! Num. of radial mesh k-points

CONTAINS

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

subroutine iofk( kf, kg, ikf, ikg )

! Finds indexes ikf and ikg such that kfmesh(ikf) <= kf < kfmesh(ikf+1)
! and kgmesh(ikg) <= kg < kgmesh(ikg+1) in logarithmic meshes of the form
! k(ik) = b*(exp((ik-1)*a)-1)

  implicit none
  real(dp), intent(in) :: kf, kg
  integer,  intent(out):: ikf, ikg

  real(dp),parameter :: ktol = 1.e-12_dp  ! 'out of range' tolerance
  real(dp),parameter :: amin = 1.e-12_dp  ! tiny denominator to avoid /0
  real(dp),save:: akf, akg, bkf, bkg
  logical, save:: first_call = .true.

  if (first_call) then
    akf = log( (kfmesh(nkf)-kfmesh(nkf-1)) / (kfmesh(2)-kfmesh(1)) ) / (nkf-2)
    akg = log( (kgmesh(nkg)-kfmesh(nkg-1)) / (kgmesh(2)-kgmesh(1)) ) / (nkg-2)
    akf = max( akf, amin )
    akg = max( akg, amin )
    bkf = (kfmesh(2) - kfmesh(1)) / (exp(akf) - 1)
    bkg = (kgmesh(2) - kgmesh(1)) / (exp(akg) - 1)
    first_call = .false.
  end if

  if (kf<kfmesh(1)-ktol .or. kf>kfmesh(nkf)+ktol .or. &
      kg<kgmesh(1)-ktol .or. kg>kgmesh(nkg)+ktol) then
     call die('vv_vdw_xc ERROR: (kf,kg) out of range')
  endif

  ikf = 1 + log( 1 + (kf-kfmesh(1))/bkf ) / akf
  ikg = 1 + log( 1 + (kg-kfmesh(1))/bkg ) / akg
  ikf = max( 1, ikf )
  ikg = max( 1, ikg )
  ikf = min( nkf-1, ikf )
  ikg = min( nkg-1, ikg )

end subroutine iofk

! -----------------------------------------------------------------------------

subroutine pofk( kf0, kg0, p0, dpdkf0, dpdkg0 )

! Finds the values and derivatives, at (kf0,kg0), of the bicubic polynomials 
! p_i(kf,kg) such that
!    y(kf,kg) = Sum_i p_i(kf,kg) * y_i
! is the bicubic spline interpolation at (kf,kg) of (any) function y with
! values y_i at mesh points (kfmesh,kgmesh)_i

  implicit none
  real(dp),intent(in) :: kf0, kg0 ! point at which the polynomials are required
  real(dp),intent(out):: p0(nkfg)      ! polynomial values at (kf0,kg0)
  real(dp),intent(out):: dpdkf0(nkfg)  ! dp/dkf at (kf0,kg0)
  real(dp),intent(out):: dpdkg0(nkfg)  ! dp/dkg at (kf0,kg0)

  integer :: ikf, ikg, ikfg, ikf0, ikg0
  real(dp):: akf, akg, bkf, bkg, dkf, dkg, &
             pkf0(nkf), dpkfdkf0(nkf), pkg0(nkg), dpkgdkg0(nkg)
  logical, save :: first_call=.true.
  real(dp),save :: pkf(nkf,nkf), d2pkfdkf2(nkf,nkf)
  real(dp),save :: pkg(nkg,nkg), d2pkgdkg2(nkg,nkg)

! Set up spline polynomial basis
  if (first_call) then
    pkf = 0
    do ikf = 1,nkf
      pkf(ikf,ikf) = 1
      call generate_spline( kfmesh, pkf(:,ikf), nkf, d2pkfdkf2(:,ikf) )
!      call generate_spline( kfmesh, pkf(:,ikf), nkf, d2pkfdkf2(:,ikf), &
!                            0._dp, 0._dp )
    end do
    pkg = 0
    do ikg = 1,nkg
      pkg(ikg,ikg) = 1
      call generate_spline( kgmesh, pkg(:,ikg), nkg, d2pkgdkg2(:,ikg) )
!      call generate_spline( kgmesh, pkg(:,ikg), nkg, d2pkgdkg2(:,ikg), &
!                            0._dp, 0._dp )
    end do
    first_call = .false.
  end if

! Find interval of qmesh in which q0 is included
  if (kf0>kfmesh(nkf) .or. kg0>kgmesh(nkg)) then   ! (kf,kg) out of range
    p0 = 0
    dpdkf0 = 0
    dpdkg0 = 0
    return
  end if
  call iofk( kf0, kg0, ikf0, ikg0 )

! Evaluate pkf polynomials of spline basis at kf0
  dkf = kfmesh(ikf0+1) - kfmesh(ikf0)
  akf = (kfmesh(ikf0+1) - kf0) / dkf   ! dadkf0 = -1/dkf
  bkf = (kf0 - kfmesh(ikf0)) / dkf     ! dbdkf0 = +1/dkf
  do ikf = 1,nkf
    pkf0(ikf) = akf*pkf(ikf0,ikf) + bkf*pkf(ikf0+1,ikf) &
              + ((akf**3-akf)*d2pkfdkf2(ikf0,ikf) +     &
                 (bkf**3-bkf)*d2pkfdkf2(ikf0+1,ikf)) * dkf**2/6
    dpkfdkf0(ikf) = - (pkf(ikf0,ikf) - pkf(ikf0+1,ikf)) / dkf &
                    - ((3*akf**2-1)*d2pkfdkf2(ikf0,ikf) -     &
                       (3*bkf**2-1)*d2pkfdkf2(ikf0+1,ikf)) * dkf/6
  end do

! Evaluate pkg polynomials of spline basis at kg0
  dkg = kgmesh(ikg0+1) - kgmesh(ikg0)
  akg = (kgmesh(ikg0+1) - kg0) / dkg   ! dadkg0 = -1/dkg
  bkg = (kg0 - kgmesh(ikg0)) / dkg     ! dbdkg0 = +1/dkg
  do ikg = 1,nkg
    pkg0(ikg) = akg*pkg(ikg0,ikg) + bkg*pkg(ikg0+1,ikg) &
              + ((akg**3-akg)*d2pkgdkg2(ikg0,ikg) +     &
                 (bkg**3-bkg)*d2pkgdkg2(ikg0+1,ikg)) * dkg**2/6
    dpkgdkg0(ikg) = - (pkg(ikg0,ikg) - pkg(ikg0+1,ikg)) / dkg &
                    - ((3*akg**2-1)*d2pkgdkg2(ikg0,ikg) -     &
                       (3*bkg**2-1)*d2pkgdkg2(ikg0+1,ikg)) * dkg/6
  end do

! Evaluate pkf*pkg polynomials at (kf0,kg0)
  ikfg = 0
  do ikg = 1,nkg
    do ikf = 1,nkf
      ikfg = ikfg+1
      p0(ikfg) = pkf0(ikf)*pkg0(ikg)
      dpdkf0(ikfg) = dpkfdkf0(ikf)*pkg0(ikg)
      dpdkg0(ikfg) = pkf0(ikf)*dpkgdkg0(ikg)
    end do
  end do

end subroutine pofk

!-----------------------------------------------------------------------------

subroutine kofn( n, gn, kf, kg, dkfdn, dkfdgn, dkgdn, dkgdgn )

! Finds Fermi and gradient wavevectors from density and gradient

  implicit none
  real(dp), intent(in) :: n          ! Electron density
  real(dp), intent(in) :: gn(3)      ! Density gradient
  real(dp), intent(out):: kf         ! Local Fermi wavevector
  real(dp), intent(out):: kg         ! |grad(n)|/n
  real(dp), intent(out):: dkfdn      ! dkf/dn
  real(dp), intent(out):: dkfdgn(3)  ! dkf/dgn
  real(dp), intent(out):: dkgdn      ! dkg/dn
  real(dp), intent(out):: dkgdgn(3)  ! dkg/dgn

  real(dp):: gn2, pi

! Trap exception for zero density
  if (n <= 1.e-15_dp) then
    kf = 0
    kg = 0
    dkfdn = 0
    dkfdgn = 0
    dkgdn = 0
    dkgdgn = 0
    return
  end if

! Find kf and kg
  pi = acos(-1._dp)
  kf = (3*pi**2*n)**(1._dp/3)
  gn2 = sum(gn**2)
  kg = sqrt(gn2)/n

! Find derivatives
  dkfdn = kf/n/3
  dkfdgn = 0
  dkgdn = -kg/n
  dkgdgn = kg*gn/gn2

end subroutine kofn

!-----------------------------------------------------------------------------

subroutine saturate( x, xc, y, dydx )

  ! Defines a function y(x,xc) = 1/(1/x^nsat+1/xc^nsat)^(1/nsat), where nsat
  ! is an integer set in the module header. This function is approximately
  ! equal to x for x<xc and it saturates to xc when x->infinity

  implicit none
  real(dp),intent(in) :: x     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: y     ! Function value
  real(dp),intent(out):: dydx  ! Derivative dy/dx

  real(dp):: z

  z = 1/x**nsat - 1/xc**nsat
  y = 1/z**(1._dp/nsat)
  dydx = y/z/x**(nsat+1)

end subroutine saturate

!-----------------------------------------------------------------------------

subroutine saturate_inverse( y, xc, x, dydx )

! Finds the inverse of the function defined in saturate subroutine:
!   y=1/(1/x^n+1/xc^n)^(1/n)  =>  x=1/(1/y^n-1/xc^n)^(1/n)

  implicit none
  real(dp),intent(in) :: y     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: x     ! Inverse function value
  real(dp),intent(out):: dydx  ! Derivative dy/dx

  real(dp):: z

  if (y<0._dp .or. y>xc) stop 'vdw:saturate_inverse: y out of range'

  z = 1/y**nsat - 1/xc**nsat
  x = 1/z**(1._dp/nsat)
  dydx = y/z/x**(nsat+1)

end subroutine saturate_inverse

!-----------------------------------------------------------------------------

subroutine set_phi_table()

! Finds and stores in memory the interpolation table (mesh points and 
! function values) for the kernel phi(k1,k2,k).

  implicit none
  integer :: ik, ikf1, ikf2, ikg1, ikg2, ik1, ik2, ir
  real(dp):: dkdk0, dphidk0, dphidkmax, dphidr0, dphidrmax, &
             k, kf1, kf2, kg1, kg2, pi, r(0:nr)

! Check that table was not set yet
  if (phi_table_set) return

! Check that kf, kg, r, and k meshes have been set
  if (.not.kmesh_set) call set_kmesh()
  if (.not.kcut_set) call die('vv_set_phi_table ERROR: kcut not set')
  forall(ir=0:nr) r(ir) = ir*dr
  pi = acos(-1.0_dp)

! Loop on (k1,k2) mesh points
  do ikg2 = 1,nkg                       ! loop on kg2
    do ikf2 = 1,nkf                     ! loop on kf2
      ik2 = ikf2 + nkf*(ik2-1)      ! combined (ikf2,ikg2) index
      do ikg1 = 1,nkg                   ! loop on kg1
        do ikf1 = 1,nkf                 ! loop on kf1
          ik1 = ikf1 + nkf*(ik1-1)  ! combined (ikf1,ikg1) index

          ! Saturated (kf,kg) values
!          kf1 = kfmesh(ikf1)
!          kf2 = kfmesh(ikf2)
!          kg1 = kgmesh(ikg1)
!          kg2 = kgmesh(ikg2)

          ! Find original (unsaturated) kf anf kg values
          call saturate_inverse( kfmesh(ikf1), kfcut, kf1, dkdk0 )
          call saturate_inverse( kfmesh(ikf2), kfcut, kf2, dkdk0 )
          call saturate_inverse( kgmesh(ikg1), kgcut, kg1, dkdk0 )
          call saturate_inverse( kgmesh(ikg2), kgcut, kg2, dkdk0 )

          ! Find kernel as a function of r12
          do ir = 0,nr
            phir(ir,ik1,ik2) = vv_phi( kf1, kf2, kg1, kg2, r(ir) )
          end do

          ! Kill kernel smoothly at r=rcut
          do ir = 0,nr
            phir(ir,ik1,ik2) = phir(ir,ik1,ik2) * cutoff( r(ir)/rcut )
          end do

          ! Find kernel in reciprocal space
          call radfft( 0, nr, rcut, phir(:,ik1,ik2), phik(:,ik1,ik2) )
          phik(:,ik1,ik2) = phik(:,ik1,ik2) * (2*pi)**1.5_dp

          ! Filter out above kcut
          phik(nk:nr,ik1,ik2) = 0

          ! Soft filter below kcut
          do ik = 1,nk
            k = ik * dk
            phik(ik,ik1,ik2) = phik(ik,ik1,ik2) * cutoff(k/kcut)
          end do

          ! Find filtered kernel in real space
          call radfft( 0, nr, kmax, phik(:,ik1,ik2), phir(:,ik1,ik2) )
          phir(:,ik1,ik2) = phir(:,ik1,ik2) / (2*pi)**1.5_dp

          ! Set up spline interpolation tables
          dphidr0 = 0
          dphidrmax = 0
          dphidk0 = 0
          dphidkmax = 0
          call spline( dr, phir(:,ik1,ik2), nr+1, dphidr0, dphidrmax, &
                       d2phidr2(:,ik1,ik2) )
          call spline( dk, phik(:,ik1,ik2), nk+1, dphidk0, dphidkmax, &
                       d2phidk2(:,ik1,ik2) )

          ! Fill symmetric elements
          phir(:,ik2,ik1) = phir(:,ik1,ik2)
          phik(:,ik2,ik1) = phik(:,ik1,ik2)
          d2phidr2(:,ik2,ik1) = d2phidr2(:,ik1,ik2)
          d2phidk2(:,ik2,ik1) = d2phidk2(:,ik1,ik2)

!          if (.false. .and. ik1==ik2) then
!            print*, 'vv_vdw_set_kcut: ik1,ik2=', ik1, ik2
!            call window( 0._dp, 5._dp, -1._dp, 4._dp, 0 )
!            call axes( 0._dp, 1._dp, 0._dp, 1._dp )
!            call plot( nr+1, r, phi, phir(:,ik1,ik2) )
!            call window( 0._dp, 10._dp, -0.05_dp, 0.15_dp, 0 )
!            call axes( 0._dp, 1._dp, 0._dp, 0.05_dp )
!            call plot( nr+1, r, r**2*phi, r**2*phir(:,ik1,ik2))
!            call show()
!          end if

        end do ! ikf1
      end do ! ikg1
    end do ! ikf2
  end do ! ikg2

! Mark table as set
  phi_table_set = .true.

end subroutine set_phi_table

! -----------------------------------------------------------------------------

subroutine set_kmesh()

! Sets mesh of q values

  implicit none
  integer :: nmesh

  if (.not.kmesh_set) then
    call set_mesh( nkf, xmax=kfcut, dxndx1=dkmaxdkmin )
    call get_mesh( nkf, nmesh, kfmesh )
    call set_mesh( nkg, xmax=kgcut, dxndx1=dkmaxdkmin )
    call get_mesh( nkg, nmesh, kgmesh )
    kmesh_set = .true.
  end if

#ifdef DEBUG_XC
  write(udebug,'(/,a,/,(10f8.4))') 'vdw:set_kmesh: kfmesh =', kfmesh
  write(udebug,'(/,a,/,(10f8.4))') 'vdw:set_kmesh: kgmesh =', kgmesh
#endif /* DEBUG_XC */

end subroutine set_kmesh

! -----------------------------------------------------------------------------

real(dp) function vv_phi( kf1, kf2, kg1, kg2, r12 )

! vdW energy kernel of Vydrov-vanVoorhis, eq.(2) of JCP 133, 244103 (2010)
! In practice, vv_phi returns n1*n2*phi(k1,k2,r12), which is smoother
! Input:
!   real(dp):: kf1, kf2 ! Fermi wavevectors at points 1 and 2
!   real(dp):: kg1, kg2 ! |grad(n)|/n at points 1 and 2
!   real(dp):: r12      ! distance between points 1 and 2

! Arguments
  implicit none
  real(dp),intent(in) :: kf1, kf2, kg1, kg2, r12

! Internal parameters and variables
  real(dp):: g1, g2, kappa1, kappa2, n1, n2, pi, w01, w02, wg1, wg2, wp1, wp2

! Find kernel
  pi = acos(-1.0_dp)
  n1 = kf1**3/(3*pi**2)           ! electron density at point 1
  n2 = kf2**3/(3*pi**2)           ! electron density at point 2
  wp1 = sqrt(4*pi*n1)             ! local plasma frequency at point 1
  wp2 = sqrt(4*pi*n2)             ! local plasma frequency at point 2
  wg1 = sqrt(vv_C*kg1**4)         ! local band gap at point 1
  wg2 = sqrt(vv_C*kg2**4)         ! local band gap at point 2
  kappa1 = vv_b*kf1**2/wp1        ! local VV kappa variable (eq.(9)) at point 1
  kappa2 = vv_b*kf2**2/wp2        ! kappa variable at point 2
  w01 = sqrt(wg1**2+wp1**2/3)     ! local w0 frequency (eq.(5)) at point 1
  w02 = sqrt(wg2**2+wp2**2/3)     ! local w0 frequency at point 2
  g1 = w01*r12**2 + kappa1        ! local g variable (eq.(3)) at point 1
  g2 = w02*r12**2 + kappa2        ! local g variable at point 2
  vv_phi = -1.5_dp/g1/g2/(g1+g2)  ! VV kernel phi (eq.(2))

! Find and return whole integrand of eq.(1)
  vv_phi = n1*n2*vv_phi

end function vv_phi

! -----------------------------------------------------------------------------

real(dp) function vv_vdw_beta()

! Returns parameter beta of VV functional

  vv_vdw_beta = vv_beta

end function vv_vdw_beta

!-----------------------------------------------------------------------------

subroutine vv_vdw_get_kmesh( nkf, nkg, kf, kg )

! Returns size and values of (kf,kg) mesh

  implicit none
  integer,          intent(out) :: nkf    ! Number of kf mesh points
  integer,          intent(out) :: nkg    ! Number of kg mesh points
  real(dp),optional,intent(out) :: kf(:)  ! Values of kf mesh points
  real(dp),optional,intent(out) :: kg(:)  ! Values of kg mesh points
  integer:: nmax
  if (.not.kmesh_set) call set_kmesh()
  if (present(kf)) then
    nmax = max( nkf, size(kf) )
    kf(1:nmax) = kfmesh(1:nmax)
  end if
  if (present(kg)) then
    nmax = max( nkg, size(kg) )
    kg(1:nmax) = kgmesh(1:nmax)
  end if
end subroutine vv_vdw_get_kmesh

! -----------------------------------------------------------------------------

subroutine vv_vdw_phi( k, phi, dphidk )

! Finds by interpolation phi(k1,k2,k) (Fourier transform of phi(k1,k2,r)) 
! for all values of k1=(kf1,kg1) and k2=(kf2,kg2) in qmesh. If the 
! interpolation table does not exist, it is calculated in the first call
! to vv_vdw_phi. It requires a previous call to vv_vdw_set_kc to set k mesh.

  implicit none
  real(dp),intent(in) :: k            ! Modulus of actual k vector
  real(dp),intent(out):: phi(:,:)     ! phi(k1,k2,k) at given k
                                      ! for all k1,k2 in (kf,kg) mesh
  real(dp),intent(out):: dphidk(:,:)  ! dphi(k1,k2,k)/dk at given k

  integer :: ik, ik1, ik2
  real(dp):: a, a2, a3, b, b2, b3

  if (.not.kcut_set) stop 'vdw_phi: ERROR: kcut must be previously set'

! Check argument sizes
  if (size(phi,1)<nkfg .or. size(phi,2)<nkfg) &
    call die('vv_vdw_phi: ERROR: size(phi) too small')

! Find phi values at point k
  if (k >= kcut) then
    phi(:,:) = 0
  else
    ! Expand interpolation inline since this is the hottest point in VdW
    ik = k/dk
    a = ((ik+1)*dk-k)/dk
    b = 1 - a
    a2 = (3*a**2 -1) * dk / 6
    b2 = (3*b**2 -1) * dk / 6
    a3 = (a**3 - a) * dk**2 / 6
    b3 = (b**3 - b) * dk**2 / 6
    do ik2 = 1,nkfg
      do ik1 = 1,ik2
!        call splint( dk, phik(:,ik1,ik2), d2phidk2(:,ik1,ik2), &
!                      nk+1, k, phi(ik1,ik2), dphidk(ik1,ik2), pr )
        phi(ik1,ik2) = a*phik(ik,ik1,ik2) + b*phik(ik+1,ik1,ik2) &
                + a3*d2phidk2(ik,ik1,ik2) + b3*d2phidk2(ik+1,ik1,ik2)
        dphidk(ik1,ik2) = (-phik(ik,ik1,ik2) &
                               +phik(ik+1,ik1,ik2) )/dk &
                - a2*d2phidk2(ik,ik1,ik2) + b2*d2phidk2(ik+1,ik1,ik2)
        phi(ik2,ik1) = phi(ik1,ik2)
        dphidk(ik2,ik1) = dphidk(ik1,ik2)
      end do
    end do
  end if

end subroutine vv_vdw_phi

!-----------------------------------------------------------------------------

subroutine vv_vdw_phiofr( r, phi )

! Finds phi(k1,k2,r) with k1=(kf1,kg1), k2=(kf2,kg2) for mesh values of
! kf's (Fermi wavevectors) and kg's (grad(n)/n)

  implicit none
  real(dp),intent(in) :: r
  real(dp),intent(out):: phi(:,:)

  integer :: ikf1, ikf2, ikg1, ikg2, ik1, ik2
  real(dp):: dphidr

  if (size(phi,1)<nkfg .or. size(phi,2)<nkfg) &
    stop 'vv_phiofr: ERROR: size(phi) too small'
  if (.not.phi_table_set) call set_phi_table()

  if (r >= rcut) then
    phi(:,:) = 0
  else
    do ik2 = 1,nkfg
      do ik1 = 1,ik2
        call splint( dr, phir(:,ik1,ik2), d2phidr2(:,ik1,ik2), &
                     nr+1, r, phi(ik1,ik2), dphidr )
        phi(ik2,ik1) = phi(ik1,ik2)
      end do ! ik1
    end do ! ik2
  end if ! (r>=rcut)

end subroutine vv_vdw_phiofr

!-----------------------------------------------------------------------------

subroutine vv_vdw_set_kcut( kc )

! Sets the reciprocal planewave cutoff kc of the integration grid, and finds 
! the interpolation table to be used by vdw_phi to obtain the vdW kernel phi
! at the reciprocal mesh wavevectors.

  implicit none
  real(dp),intent(in):: kc  ! Planewave cutoff: k>kcut => phi=0

  real(dp):: pi

  ! Set kcut and radial mesh parameters, if not already set
  if (.not. kc==kcut) then
    pi = acos(-1._dp)
    dr = rcut / nr
    dk = pi / rcut
    kmax = pi / dr
    nk = int(kc/dk) + 1
    if (nk>nr) stop 'vv_vdw_set_kcut: ERROR: nk>nr'
    kcut = kc
    kcut_set = .true.
  end if

  ! Set (kf,kg) mesh and phi table
  call set_kmesh()
  call set_phi_table()

#ifdef DEBUG_XC
  write(udebug,'(a,5f8.3)') 'vv_vdw_set_kcut: kfcut,kgcut,rcut,kcut,kmax=', &
    kfcut, kgcut, rcut, kc, kmax
#endif /* DEBUG_XC */

end subroutine vv_vdw_set_kcut

!-----------------------------------------------------------------------------

subroutine vv_vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )

! Finds the value and derivatives of theta_i(rho,grad_rho) of eq.(8)
! of Roman-Soler PRL 2009. In practice, theta_i=p_i rather than rho*p_i,
! beacuse p_i is used here to expand rho1*rho2*phi rather than only phi.

  implicit none
  integer, intent(in) :: nspin                 ! Number of spin components
  real(dp),intent(in) :: rhos(nspin)           ! Electron spin density
  real(dp),intent(in) :: grhos(3,nspin)        ! Spin density gradient
  real(dp),intent(out):: theta(nkfg)           ! Exp. polynomials at (kf,kg)
  real(dp),intent(out):: dtdrho(nkfg,nspin)    ! dtheta(iq)/drhos
  real(dp),intent(out):: dtdgrho(3,nkfg,nspin) ! dtheta(iq)/dgrhos(ix)

  integer :: ikfg, is, ns
  real(dp):: rho, grho(3), dpdkf(nkfg), dpdkg(nkfg), dkfdrho, dkfdgrho(3), &
             dkgdrho, dkgdgrho(3), kf, kg, p(nkfg)

  ! Sum spin components of electron density
  ns = min(nspin,2)     ! num. of diagonal spin components (if noncollinear)
  rho = sum(rhos(1:ns))             ! local electron density
  grho = sum(grhos(:,1:ns),dim=2)  ! local density gradient

  ! Find local Fermi and gradient wavevectors, and their derivatives
  call kofn( rho, grho, kf, kg, dkfdrho, dkfdgrho, dkgdrho, dkgdgrho )

  ! Find expansion polynomials of integrand kernel of nonlocal vdW energy
  call pofk( kf, kg, p, dpdkf, dpdkg )

  ! Find theta functions and their derivatives with respect to rho and grho
  do ikfg = 1,nkfg
    theta(ikfg) = p(ikfg)  ! because we expand rho1*rho2*phi rather than phi
    do is = 1,ns
      dtdrho(ikfg,is) = dpdkf(ikfg)*dkfdrho + dpdkg(ikfg)*dkgdrho
      dtdgrho(1:3,ikfg,is) = dpdkf(ikfg)*dkfdgrho(1:3) &
                           + dpdkg(ikfg)*dkgdgrho(1:3)
    end do
  end do

end subroutine vv_vdw_theta

END MODULE m_vv_vdwxc
