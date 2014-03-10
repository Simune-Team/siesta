!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp9in
!! NAME
!! psp9in
!!
!! FUNCTION
!! Initialize pspcod=9 (Pseudopotentials from the XML format):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG, AF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=maximum number of points in real space grid in the psp file
!!   angular momentum of nonlocal pseudopotential
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optnlxccc=option for nl XC core correction (input variable)
!!  qgrid(mqgrid)=values of q (or |G|) on grid from 0 to qmax
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r} dr]$ (hartree)
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  qchrg is not used, and could be suppressed later
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr)
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!
!! PARENTS
!!      pspatm_abinit
!!
!! CHILDREN
!!      atomxc,cc_derivatives,close_xml_t,open_xml_file,parse,psp5lo,psp5nl
!!      schro_eq,smoothvlocal,spline,splint,vhrtre,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psp9in(filpsp,ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mqgrid,nproj,n1xccc,qchrg,qgrid,&
&                  useylm,vlspl,xcccrc,xccc1d,zion,znucl)

 use m_profiling

 use defs_basis
 use m_splines
 use m_schrodinger,    only: schro_eq
#if defined HAVE_TRIO_FOX
 use fox_sax
 use m_xml_pseudo_types
 use m_xml_pseudo
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp9in'
 use interfaces_14_hidewrite
 use interfaces_21_psiesta_noabirule
 use interfaces_65_psp, except_this_one => psp9in
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mpsang,mqgrid,n1xccc
 integer,intent(in) :: useylm
 integer,intent(out) :: mmax
 real(dp),intent(in) :: zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc
 character(len=fnlen),intent(in) :: filpsp
!arrays
 integer,intent(out) :: indlmn(6,lmnmax),nproj(mpsang)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: ekb(lnmax),ffspl(mqgrid,2,lnmax),vlspl(mqgrid,2)
 real(dp),intent(out) :: xccc1d(n1xccc,6)

!Local variables-------------------------------
!no_abirules
#if defined HAVE_TRIO_FOX
 integer, parameter :: nrmax = 20000   ! Maximum number of points 
                                       !   in the functions read from the 
                                       !   XML pseudopotential file
 integer, parameter :: nkbmx = 2       ! Maximum number of Kleinman-Bylander 
                                       !   projectors for each angular momentum
 real(dp), parameter :: Rmax_kb_default = 60.0_dp ! Radius of the sphere where
                                                  !   the atomic wavefunctions
                                                  !   will be computed
 character*10  :: functl   
 character*10  :: author  
 integer :: ipsang,ir,jj,jpsang,maxn_pots,mm,npts
 integer :: il,index,irelt,nrwf,nnodes,nprin,l,ii
 real(dp) :: ainc,al,amesh,cmax,dnorm,fchrg,phi,r2,ratio,dc,zval
 real(dp) :: rmax,rmin,scale,step,yp1,ypn,z,chgvps,rpb,ea,rmaxwf,dnrm
 real(dp) :: ex, ec, dx, rgauss, rgauss2
 logical :: vlocsmooth,nicore,use_read_wfnl
 character(len=500) :: message
 real(dp),allocatable :: ekb_tmp(:),ffspl_tmp(:,:,:),funrad(:),rad(:)
 real(dp),allocatable :: ratm(:),rfhi(:),vloc(:),vps(:,:),vpspll(:,:)
 real(dp),allocatable :: drdi(:),metric(:)
 real(dp),allocatable :: ve(:),vxc(:),rho(:),auxrho(:),chcore(:)
 real(dp),allocatable ::       ff1(:),ff2(:)
 real(dp),allocatable :: vlocal(:),wfll(:,:),work_spl(:)
 real(dp),allocatable :: wfllatm(:,:), erefkb(:,:)
! integer :: ii,irad,i1xccc,lhigh,ll,mmax2,pspcod
!real(dp) :: al_announced,der1,dern,factor,ff1fact
!real(dp),allocatable ::       ff3(:),ff4(:)
!real(dp),allocatable :: gg(:),gg1(:),gg2(:),gg3(:),gg4(:),radbis(:)
!----------------------------------------------------------------------------
! character*10 functl  : Exchange and correlation functional
! character*10 author  : Exchange and correlation parametrization
! real(dp) z    : Atomic number of the element.
! real(dp) rmin : Value to determine the innermost point of the mesh
!                 in the FHI code ( rfhi(1) = rmin/dble(atomic-number) )
! real(dp) rmax : Outermost point of the mesh in the FHI code.
! real(dp) ainc : Default value for the ratio for the geometrical sequence
!                 used to define the grid in the FHI code.
! real(dp) amesh: Actual value for the ratio for the geometrical sequence
!                 used to define the grid in the FHI code.
! integer  mmax : Maximum number of points in the radial grid.
! real(dp) rfhi : Points of the radial grid in the FHI code,
!                 rfhi(ir) = rfhi(ir-1) * amesh .
! real(dp) scale: Scale used to define the points of the logarithmic
!                 radial grid of the ATM3 code.
! real(dp) step : Step used to define the points of the logarithmic
!                 radial grid of the ATM3 code.
! integer  npts : Number of points in the logarithmic radial grid of
!                 the ATM3 code.
! real(dp) ratm : Points of the radial grid in the ATM code,
!                 ratm(ir) = scale *  [ exp( step*ir) - 1 ]
! real(dp) drdi : Derivative of the radial distance respect the mesh index
! real(dp) metric : Metric array
! real(dp) rad    : Radial grid (geometrical sequence).
! real(dp) vpspll : Semilocal component of the pseudopotential.
!                   Units in Hartress.
!                   First argument: index in the radial grid.
!                   Second argument: angular momentum.
! real(dp) wfll   : Radial part of the pseudowave functions
!                   directly read from the pseudopotential file.
!                   wfll = u_n,l (r) = 1/r R_n,l(r)
!                   Units: electrons/bohr**(3/2).
!                   Normalized to 1.
!                   First argument: index in the radial grid.
!                   Second argument: angular momentum.
! real(dp) wfllatm: Radial part of the pseudowave functions, after
!                   solving the Schrodinger equation with the semilocal
!                   components of the pseudopotential.
!                   wfllatm = u_n,l (r) = 1/r R_n,l(r)
!                   Units: electrons/bohr**(3/2).
!                   Normalized to 1.
!                   First argument: index in the radial grid.
!                   Second argument: angular momentum.
! real(dp) erefkb : Reference energy (in Ry) for the calculation of the KB proj
! real(dp) chgvps : Valence charge of the pseudoion for which the pseudo
!                   was generated in the ATM code
!                   (it might not correspond with the nominal valence charge
!                   of the atom if the pseudo has been generated for an ionic
!                   configuration, for instance when semicore has been 
!                   explicitly included in the valence).
!                   For instance, for Ba with the semicore in valence,
!                   (atomic reference configuration 5s2 5p6 5d0 4f0),
!                   chgvps = 8  (two in the 5s and six in the 5p)
!                   zval   = 10 (two in the 5s, six in the 5p, and two in the
!                                6s. These two last not included in the 
!                                reference atomic configuration)
! real(dp) rho    : Valence charge density generated by the ATM code.
!                   As read from the ATM file, it is angularly integrated
!                   (i.e. multiplied by 4*pi*r^2).
! real(dp) ve     : Electrostatic potential generated by the valence charge
!                   density, readed from the pseudo file
! real(dp) vxc    : Exchange and correlation potential
!                   density, readed from the pseudo file
! real(dp) chcore : Core charge density generated by the ATM code.
!                   As read from the ATM file, it is angularly integrated
!                   (i.e. multiplied by 4*pi*r^2).
! real(dp) auxrho : Sum of the valence charge and core charge (if NLCC included)
!                   densities to compute the atomic exchange and correlation
!                   potential. 
!                   auxrho is NOT angularly integrated 
!                   (i.e. not multiplied by 4*pi*r^2).
! integer irelt   : Flag that determines whether the atomic calculation to 
!                   generate the pseudopotential was relativistic (irelt = 1)
!                   or no relativistic (irelt = 0)
! logical nicore  : Flag that determines whether non-linear-core-corrections
!                   are considered or not
! integer nrwf    : Number of points in the ATM logarithmic grid required
!                   to store the pseudowavefunctions
! logical use_read_wfnl : Determine whether the code will use the 
!                         wave functions generated within the subroutine or
!                         those directly read from the XML file.
! logical vlocsmooth : Flag to determine whether one of the semilocal components
!                      of the pseudopotential is used as local part, or 
!                      if the Siestalike smooth local component is generated
!                      To facilitate the comparison with Siesta, it is set to 
!                      true
! real(dp) rgauss : Radius at which the semilocal components of the 
!                   pseudopotentials have all converged to the same 
!                   tail as the last one
! real(dp) rgauss2: Radius at which all the pseudopotentials 
!                   have converged to 2*z/r
!---------------------------------------------------------------------------
! Variables for the non-linear core corrections ------------------------
!integer              :: nptscore
!----------------------------------------------------------------------------
! xcccrc           : XC core correction cutoff radius (bohr)
!                    It is defined as the radius where the pseudo-core
!                    charge density becomes zero
!                    (here we hae set up a tolerance of 1.d-12).
!----------------------------------------------------------------------------
 integer                  :: iostat
 type(xml_t)              :: fxml
 type(pseudo_t), pointer  :: psxml
#endif

! ***************************************************************************

 if(.false.)write(std_out,*)filpsp ! Just to keep filpsp when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)lloc   ! Just to keep lloc when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)lmax   ! Just to keep lmax when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)qgrid  ! Just to keep qgrid when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)useylm ! Just to keep useylm when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)zion   ! Just to keep zion when HAVE_TRIO_FOX is false
 if(.false.)write(std_out,*)znucl  ! Just to keep znucl when HAVE_TRIO_FOX is false
#if defined HAVE_TRIO_FOX
!DEBUG
!write(std_out,*)' psp9in : enter and stop '
!stop
!ENDDEBUG

!Determine whether the code will use the wave functions previously
!generated or those directly read from the XML file.
!For a matter of compatibility with the wave functions used in 
!Siesta to generate the KB projectors, 
!we will use by default the wave functions generated within this
!subroutine.
!But it can be trivially changed setting the following variable to .true.
 use_read_wfnl = .false.

!Determine whether the smooth local part of the pseudopotential
!generated "a-la-Siesta" will be used
!see the Siesta technical paper, 
!J. M. Soler et al. J. Phys.: Condens. Matter 14, 2745 (2002)
!If this variable is set to .false. a semilocal component of the 
!pseudopotential us used as local part
 vlocsmooth = .true.

 call open_xml_file(fxml,filpsp,iostat)
 if (iostat /=0) stop "Cannot open file"

 call parse(fxml,pcdata_chunk,startElement_handler=begin_element,endElement_handler=end_element)

 call close_xml_t(fxml)
 if (iostat /=0) stop "Cannot close file"

 psxml => pseudo

!Define the mesh used by the FHI pseudopotential code
!The atomic number of the element is read from the header of the XML file
!The default values for rmin, rmax, and ainc are:
!These are taken from the "fine mesh usually overdoing somewhat"
!of the file default.h of the FHI98PP code
 rmin = 0.00625_dp
 rmax = 80.0_dp
 ainc = 1.0123_dp
 z    = psxml%header%atomicnumber
 zval = psxml%header%zval
!---

!Determine the maximum number of points in the grid ---
 amesh    = dble(ainc)
 al       = log(amesh)
 cmax     = dble(rmax)/rmin
 mmax     = log(cmax*z)/al
 if(mod(mmax,2) .eq. 0) mmax = mmax + 1
!DEBUG
 write(std_out,*)' psp9in : parameters to define the points of the grid '
 write(std_out,*)' psp9in : rmin = ', rmin
 write(std_out,*)' psp9in : rmax = ', rmax
 write(std_out,*)' psp9in : ainc = ', ainc
 write(std_out,*)' psp9in : mmax = ', mmax
!ENDDEBUG
!---

!Allocate the radial grid in the FHI code, and define the points ---
 ABI_ALLOCATE( rfhi,(mmax))
 ABI_ALLOCATE( rad,(mmax))
 rfhi(1) = dble(rmin)/z
 do ir = 2, mmax
   rfhi(ir) = amesh*rfhi(ir-1)
 end do
 rad(:) = rfhi(:)
!! DEBUG
! do ir = 2, mmax
!   write(std_out,'(i5,f20.12)')ir, rfhi(ir)
! end do
!! ENDDEBUG
!---

!Now, define the mesh used by ATM3 pseudopotential code.
!We will assume at the present stage, that all the radial grids
!for the different magnitudes are the same
 npts  = psxml%pot(1)%V%grid%npts + 1
 scale = psxml%pot(1)%V%grid%scale
 step  = psxml%pot(1)%V%grid%step

!Allocate the radial grid in the FHI code, and define the points ---
 ABI_ALLOCATE( ratm,(npts))
 ABI_ALLOCATE( drdi,(npts))
 ABI_ALLOCATE( metric,(npts))
 do ir = 1, npts
   ratm(ir) = scale * (exp(step*(ir-1))-1)
 end do
!Calculate drdi and s
!drdi is the derivative of the radial distance respect to the mesh index
!i.e. ratm(ir)= b*[ exp( a*(i-1) ) - 1 ] and therefore
!drdi=dr/di =a*b*exp(a*(i-1))= a*[ratm(ir)+b]
 rpb = scale
 ea  = dexp(step)
 do ir = 1, npts
   drdi(ir) = step * rpb
   metric(ir)    = dsqrt( step * rpb )
   rpb      = rpb * ea
 end do
!---

!Allocate the radial variables (semilocal pseudopotentials and
!pseudoatomic orbitals)
 ABI_ALLOCATE( vpspll,(mmax,mpsang))
 ABI_ALLOCATE( wfll,(mmax,mpsang))
 ABI_ALLOCATE( vps,(npts,0:lmax))

!Interpolate atomic pseudopotential for each l, filling up arrays vpspll
!Note: put each l into vpspll(:,l+1)
 do il = 1, mpsang

   ABI_ALLOCATE( funrad,(npts))
   ABI_ALLOCATE( ff2,(npts))

   funrad(2:npts) = psxml%pot(il)%V%data(1:npts-1)
   r2        = ratm(2) / ( ratm(3) - ratm(2) )
   funrad(1) = funrad(2) - r2 * ( funrad(3) - funrad(2) )
   yp1       = ( funrad(2) - funrad(1) ) / &
&   ( ratm(2) - ratm(1) )
   ypn       = ( funrad(npts) - funrad(npts-1) )/&
&   ( ratm(npts) - ratm(npts-1) )

!  Store the semilocal components of the pseudopotential in the variable vps,
!  required to compute the local part of the pseudopotential.
!  vps is read in r*V format. Here we transform it to V format
   vps(2:npts,il-1) = funrad(2:npts) / ratm(2:npts)
   vps(1,il-1)      = vps(2,il-1)

!  !DEBUG
!  write(std_out,*)
!  write(std_out,*)'# il = ', il-1
!  do ir = 2, npts
!  write(std_out,'(3f20.12)')ratm(ir), &
!  &          psxml%pot(il)%V%data(ir-1)/(ratm(ir)*2.0_dp), funrad(ir)
!  end do
!  !ENDDEBUG

!  Be careful, the interpolation is made with funrad,
!  and not with vps.
!  funrad is still in r*V format. That is why we multiply
!  afterwards by (1/r) again.

   call spline ( ratm, funrad, npts, yp1, ypn, ff2 )

   call splint( npts, ratm, funrad, ff2, mmax, rad, vpspll(:,il) )

!  Transform the pseudopotential from r*V format to V format,
!  and from Rydbergs to Hartress

   do ir = 1, mmax
     vpspll(ir,il) = vpspll(ir,il) / ( 2.0_dp * rad(ir) )
   end do

   ABI_DEALLOCATE( funrad)
   ABI_DEALLOCATE( ff2)
 end do

!!DEBUG
!do il = 1, mpsang
!write(std_out,*)
!write(std_out,*)'# il = ', il-1
!do ir = 1, mmax
!write(std_out,'(2f20.12)')rad(ir), vpspll(ir,il)
!end do
!end do
!!ENDDEBUG

! Solve the Schrodinger equation for the isolated pseudoatom for each 
! angular momentum shell

 ABI_ALLOCATE( ve,(nrmax)     )
 ABI_ALLOCATE( vxc,(nrmax)    )
 ABI_ALLOCATE( rho,(nrmax)    )
 ABI_ALLOCATE( auxrho,(nrmax) )
 ABI_ALLOCATE( chcore,(nrmax) )

! chgvps = psxml%gen_zval

! Read the valence charge density from the ATM file.
 rho(2:npts) = psxml%valence_charge%data(1:npts-1)
 r2     = ratm(2) / ( ratm(3) - ratm(2) )
 rho(1) = rho(2) - r2 * ( rho(3) - rho(2) )

! Find the Hartree potential created by a radial electron density
! using the Numerov's method to integrate the radial Poisson equation.
! The input charge density at this point has to be angularly integrated.
 call vhrtre( rho, ve, ratm, drdi, metric, npts, step )

! Determine whether non-linear-core-corrections are included.
! If they are, read the core charge density
 select case(psxml%header%core_corrections)
   case("yes")
     nicore = .true.
!    Read the core charge density from the pseudo file
     chcore(2:npts) = psxml%core_charge%data(1:npts-1)
     r2        = ratm(2) / ( ratm(3) - ratm(2) )
     chcore(1) = chcore(2) - r2 * ( chcore(3) - chcore(2) )
   case("no")
     nicore    = .false.
     chcore(:) = 0.0_dp
 end select

! Compute the exchange and correlation potential in the atom
! Note that atomxc expects true rho(r), not 4 * pi * r^2 * rho(r)
! We use auxrho for that

 do ir = 2, npts
   r2 = ratm(ir)**2
   r2 = 4.0_dp * pi * r2
   dc = rho(ir) / r2
   if( nicore )  dc = dc + chcore(ir) / r2
   auxrho(ir) = dc
 end do
 r2        = ratm(2) / (ratm(3)-ratm(2))
 auxrho(1) = auxrho(2) - ( auxrho(3) - auxrho(2) ) * r2

! Determine whether the atomic calculation to generate the pseudopotential
! is relativistic or not
 select case(psxml%header%relativistic)
   case(.true.)
     irelt = 1
   case(.false.)
     irelt = 0
 end select 
!DEBUG
 write(std_out,*)' psp9in : relativistic pseudopotential? (1=yes, 0 =no) '
 write(std_out,*)' psp9in : irelt = ', irelt
!ENDDEBUG

! Determine the exchange and correlation functional and parametrization
 functl = psxml%header%xcfunctionaltype(1:10)
 select case(psxml%header%xcfunctionalparametrization)
   case("Ceperley-Alder")
     author = "PZ"
   case('Perdew-Burke-Ernzerhof')
     author = "PBE"
   case('Becke-Lee-Yang-Parr')
     author = "BLYP"
   case('Wu-Cohen')
     author = "WC"
 end select 

!DEBUG
 write(std_out,*)' psp9in : exchange and correlation functional ', functl 
 write(std_out,*)' psp9in : exchange and correlation parametrization ', author 
!ENDDEBUG

 call atomxc( functl, author, irelt, npts, nrmax, ratm,     &
& 1, auxrho, ex, ec, dx, dc, vxc )

! Add the exchange and correlation potential to the Hartree potential
 ve(1:npts) = ve(1:npts) + vxc(1:npts)

!!DEBUG
!do ir = 1, npts
!write(std_out,'(3f20.12)')ratm(ir), auxrho(ir), ve(ir)
!end do
!!ENDDEBUG

!
! Redefine the array s for the Schrodinger equation integration
!
 metric(2:npts) = drdi(2:npts) * drdi(2:npts)
 metric(1) = metric(2)

! Define the radius of the sphere where the atomic wavefunctions
! will be computed, and the corresponding index of the point
! in the logarithmic grid (nrwf)
 rmaxwf = Rmax_kb_default
! DEBUG        
 write(std_out,*)' psp9in : radius of the sphere where '
 write(std_out,*)' psp9in : the atomic wavefunctions will be computed '
 write(std_out,*)' psp9in : Rmax_kb = ', rmaxwf
! ENDF DEBUG   
 nrwf   = nint(log(rmaxwf/scale+1.0_dp)/step) + 1
 nrwf   = min(nrwf,npts)

 ABI_ALLOCATE( wfllatm,(nrmax,mpsang))
 ABI_ALLOCATE( erefkb,(nkbmx,0:lmax))
 erefkb(:,:) = huge(1.0_dp)
! Loop over all the angular momenta for which we are going to compute
! a KB projector
 do il = 1, mpsang
   l = il -1
!  Initialize the atomic wavefunction
   wfllatm(:,il) = 0.0_dp
   nnodes = 1
   nprin  = l+1
!  Solve the Schrodinger equation for the isolated atom.
!  The boundary condition is that the eigenstates have to be zero at r(nrwf) 
   call schro_eq( zval, ratm, vps(1,l), ve, metric, drdi,    &
&   nrwf, l, step, scale, nnodes, nprin,  &
&   erefkb(1,l), wfllatm(1,il) )
!  Normalization of the eigenstates inside a sphere of radius Rmax=nrwf
   dnrm = 0.0d0
   do ir = 1, nrwf
     dnrm =dnrm + drdi(ir) * wfllatm(ir,il)**2
   end do
   dnrm = sqrt(dnrm)
   do ir = 1, nrwf
     wfllatm(ir,il) = wfllatm(ir,il)/dnrm
   end do
!!DEBUG
!   write(std_out,'(a,i5)')' # il      = ', il
!   write(std_out,'(a,i5)')' # l       = ', l
!   write(std_out,'(a,i5)')' # nnodes  = ', nnodes
!   write(std_out,'(a,i5)')' # nprin   = ', nprin
!   do ir = 1, npts
!     write(std_out,'(3f20.12)')ratm(ir), wfllatm(ir,il)
!   end do
!   write(std_out,'(a,i5)')
!!ENDDEBUG
 end do

 ABI_DEALLOCATE( ve      )
 ABI_DEALLOCATE( vxc     )
 ABI_DEALLOCATE( rho     )
 ABI_DEALLOCATE( auxrho  )
 ABI_DEALLOCATE( erefkb  )
 ABI_DEALLOCATE( drdi     )
 ABI_DEALLOCATE( metric   )


!Interpolate the radial wave functions ---
 do il = 1, mpsang

   ABI_ALLOCATE( funrad,(npts))
   ABI_ALLOCATE( ff2,(npts))

   if( use_read_wfnl ) then
     funrad(2:npts) = psxml%pswf(il)%V%data(1:npts-1)
     r2        = ratm(2) / ( ratm(3) - ratm(2) )
     funrad(1) = funrad(2) - r2 * ( funrad(3) - funrad(2) )
   else
     funrad(1:npts) = wfllatm(1:npts,il)
   end if

   yp1       = ( funrad(2) - funrad(1) ) / &
&   ( ratm(2) - ratm(1) )
   ypn       = ( funrad(npts) - funrad(npts-1) )/&
&   ( ratm(npts) - ratm(npts-1) )

!! DEBUG
!  write(std_out,*)
!  write(std_out,*)'# il = ', il-1
!  do ir = 2, npts
!  write(std_out,'(3f20.12)')ratm(ir), psxml%pswf(il)%V%data(ir-1), funrad(ir)
!  end do
!! ENDDEBUG

   call spline ( ratm, funrad, npts, yp1, ypn, ff2 )

   call splint( npts, ratm, funrad, ff2, mmax, rad, wfll(:,il) )

   ABI_DEALLOCATE( funrad )
   ABI_DEALLOCATE( ff2 )

!  Normalize the pseudowave functions to 1.
   dnorm = 0.0_dp
   do ir = 2, mmax
     phi   = wfll(ir,il)
     dnorm = dnorm + phi**2 * ( rad(ir) - rad(ir-1) )
   end do
   dnorm = sqrt(dnorm)

   do ir = 1, mmax
     wfll(ir,il) = wfll(ir,il)/dnorm
   end do

 end do

!!DEBUG
!do il = 1, mpsang
!write(std_out,*)
!write(std_out,*)'# il = ', il-1
!do ir = 1, mmax
!write(std_out,'(2f20.12)')rad(ir), wfll(ir,il)
!end do
!end do
!!ENDDEBUG

!Take care of the non-linear core corrections

 select case(psxml%header%core_corrections)
   case("yes")
!    In Abinit, at least for the Troullier-Martins pseudopotential,
!    the pseudocore charge density and its derivatives (xccc1d)
!    are introduced in a linear grid.
!    This grid is normalized, so the radial coordinates run between
!    from 0 and 1 (from 0 to xcccrc, where xcccrc is the radius
!    where the pseudo-core becomes zero).

!    Allocate some of the arrays
     ABI_ALLOCATE( funrad,(npts))
     ABI_ALLOCATE( ff1,(npts))
     ABI_ALLOCATE( ff2,(npts))

!    Store the value of the pseudo-core charge.
!    It is read in a logarithmic grid.
     funrad(2:npts) = psxml%core_charge%data(1:npts-1)
     funrad(2:npts) = funrad(2:npts)/(4.0_dp*pi*ratm(2:npts)**2)
     r2        = ratm(2) / ( ratm(3) - ratm(2) )
     funrad(1) = funrad(2) - r2 * ( funrad(3) - funrad(2) )
     yp1       = ( funrad(2) - funrad(1) ) / &
&     ( ratm(2) - ratm(1) )
!     yp1       = zero
     ypn       = ( funrad(npts) - funrad(npts-1) )/&
&     ( ratm(npts) - ratm(npts-1) )

!    determine xcccrc where the pseudocore becomes 0
     do jj=npts,1,-1
       if (funrad(jj) > 1.0d-12) then
         xcccrc=ratm(jj)
         exit
       end if
     end do

!    Find the first derivative of the pseudo-core charge
!    in the logarithmic grid.
     ff1(1)=yp1
     do jj=2,npts-1
       ff1(jj)=(funrad(jj+1)-funrad(jj))/(ratm(jj+1)-ratm(jj))
     end do
     ff1(npts)=zero ! Could try to set this to 0 as well

!    Find the second derivative of the pseudo-core charge
!    in the logarithmic grid.
!    Be careful, this is very noisy at r->0
     call spline( ratm, funrad, npts, yp1, ypn, ff2 )

!    call cc_derivatives to get 3rd 4th and 5th derivatives,
!    and everything splined onto regular grid [0:xcccrc]
!    in xccc1d
     write(std_out,*) 'psp9in: about to call cc_derivatives'
     write(std_out,*) npts,n1xccc,xcccrc

     call cc_derivatives(ratm,funrad,ff1,ff2,npts,n1xccc,xcccrc,xccc1d)

!!    DEBUG
!    do ii = 2, npts
!    write(std_out,'(4f20.12)') &
! &    ratm(ii),chcore(ii)/(4.0_dp*pi*ratm(ii)**2),funrad(ii),ff2(ii)
!    enddo
!    
!    write(std_out,*) '# psp9in NLCC data ', n1xccc
!    do ii = 1, n1xccc
!    write(std_out,'(7e20.8)')xcccrc*(ii-1.d0)/(n1xccc-1.d0),xccc1d(ii,1),xccc1d(ii,2),xccc1d(ii,3),xccc1d(ii,4),xccc1d(ii,5),xccc1d(ii,6)
!    enddo
!!    ENDDEBUG

     ABI_DEALLOCATE( funrad )
     ABI_DEALLOCATE( ff1 )
     ABI_DEALLOCATE( ff2 )
   case("no")
     write(message, '(a)' ) '  No XC core correction.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     xcccrc = 0.0_dp ; fchrg = 0.0_dp ; qchrg = 0.0_dp
 end select

!Deallocate the different parts of psxml
 maxn_pots=size(psxml%pswf)

 do il = 1, maxn_pots
   if(associated( psxml%pswf(il)%V%data)) then
     ABI_DEALLOCATE (psxml%pswf(il)%V%data)
   end if
   if(associated( psxml%pot(il)%V%data))   then
     ABI_DEALLOCATE (psxml%pot(il)%V%data)
   end if
 end do
 if(associated( psxml%core_charge%data))   then
   ABI_DEALLOCATE (psxml%core_charge%data)
 end if
 if(associated (psxml%valence_charge%data))  then
   ABI_DEALLOCATE (psxml%valence_charge%data)
 end if

!Compute in real(dp) al : the announced amesh is inaccurate.
 ratio=rad(mmax)/rad(1)
 al=log(ratio)/dble(mmax-1)

!DEBUG
!write(std_out,*)' psp9in : al ; al_announced =',al,al_announced
!allocate(radbis(mmax))
!write(std_out,*)' psp9in : lloc  ',lloc
!do ipsang=1,lmax+1
!write(std_out,*)' psp9in : ipsang  ',ipsang
!do irad=1,mmax
!write(std_out,*)irad,rad(irad),wfll(irad,ipsang),vpspll(irad,ipsang)
!end do
!end do
!deallocate(radbis)
!ENDDEBUG

!Define the local component of the pseudopotential
 ABI_ALLOCATE(vloc,(mmax))

 ABI_ALLOCATE( vlocal,(npts))

 if ( vlocsmooth ) then
   call smoothvlocal( lmax,npts,scale,step,vlocal,vps,zval,rgauss,rgauss2 )
! DEBUG
   write(std_out,*)' psp9in : Radius at which the semilocal components of the '
   write(std_out,*)' psp9in : pseudopotentials have all converged to the same '
   write(std_out,*)' psp9in : tail as the last one '
   write(std_out,*)' psp9in : rgauss = ', rgauss
   write(std_out,*)' psp9in : Radius at which all the pseudopotentials '
   write(std_out,*)' psp9in : have converged to 2*z/r '
   write(std_out,*)' psp9in : rgauss2 = ', rgauss2
! ENDDEBUG
!  vlocal is the soft local pseudopotential used in Siesta.
!  vlocal is computed in Ry. We transform it into Hartrees.
   vlocal(1:npts) = vlocal(1:npts) / 2.0_dp
   yp1  = ( vlocal(2) - vlocal(1) ) / &
&   ( ratm(2) - ratm(1) )
   ypn  = ( vlocal(npts) - vlocal(npts-1) )/&
&   ( ratm(npts) - ratm(npts-1) )

   ABI_ALLOCATE( ff2,(npts))

   call spline( ratm, vlocal, npts, yp1, ypn, ff2 )

   call splint( npts, ratm, vlocal, ff2, mmax, rad, vloc )

! We make the local components of the pseudopotential equal to one of them
! beyond rgauss2. This ensures that (vpspll - vloc ) = 0 beyond that
! point avoinding numerical noise.

   do ir = 1, mmax
     if( rad(ir) .gt. rgauss2) vloc(ir) = vpspll(ir,1)
   end do 

!! DEBUG
!  do ir = 1, mmax
!  write(std_out,'(2f20.12)') rad(ir), vloc(ir)
!  end do
!! ENDDEBUG     

   ABI_DEALLOCATE( ff2)
 else
!  vloc(:)=Vlocal(r), lloc=0, 1, or 2 or -1 for avg.
!  Copy appropriate nonlocal psp for use as local one
   vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )
 end if

!--------------------------------------------------------------------
!Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform)
!to get q^2 V(q).

 call psp5lo(al,epsatm,mmax,mqgrid,qgrid,&
& vlspl(:,1),rad,vloc,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_spl,(mqgrid))
 call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_spl)

!--------------------------------------------------------------------
!Take care of non-local part

 ABI_ALLOCATE(ekb_tmp,(mpsang))
 ABI_ALLOCATE(ffspl_tmp,(mqgrid,2,mpsang))

!Zero out all Kleinman-Bylander energies to initialize
 ekb_tmp(:)=0.0_dp

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then
   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 else

!  ----------------------------------------------------------------------
!  Compute KB form factors and fit splines

   call psp5nl(al,ekb_tmp,ffspl_tmp,lmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,&
&   vpspll,wfll)

 end if

!Define the number of projector per angular momentum shell (1 by default)
 do ipsang = 1, lmax + 1
   nproj(ipsang) = 1
 end do
 jj=0;index=0;indlmn(:,:)=0
 do ipsang=1,lmax+1
!  nproj had been set at 1, by default
   if(abs(ekb_tmp(ipsang))<tol10)then
     nproj(ipsang)=0
   end if
!  Possible values for nproj in this routine : 0 or 1.
   if(nproj(ipsang)==1)then
     if (useylm==1) then
       jj=jj+1
       do mm=1,2*ipsang-1
         index=index+1
         indlmn(1,index)=ipsang-1
         indlmn(2,index)=mm-ipsang
         indlmn(3,index)=1
         indlmn(4,index)=mm+(ipsang-1)*(ipsang-1)
         indlmn(5,index)=jj
         indlmn(6,index)=1
       end do
     else
       jj=jj+1
       index=index+1
       indlmn(1,index)=ipsang-1
       indlmn(2,index)=0
       indlmn(3,index)=1
       indlmn(4,index)=ipsang+(ipsang-1)*(ipsang-1)
       indlmn(5,index)=jj
       indlmn(6,index)=1
     end if
   end if
 end do

!Transfer ekb and ffspl to their definitive location
 jpsang=1
 do ipsang=1,lmax+1
   if(nproj(ipsang)/=0)then
     ekb(jpsang)=ekb_tmp(ipsang)
     ffspl(:,:,jpsang)=ffspl_tmp(:,:,ipsang)
     jpsang=jpsang+1
     if(jpsang>lnmax)then
       write(message,'(6a,2i6)') ch10,&
&       ' psp9in : BUG -',ch10,&
&       '  Problem with the dimension of the ekb and ffspl arrays.',ch10,&
&       '  ipsang,lnmax=',ipsang,lnmax
     end if
   end if
 end do

 ABI_DEALLOCATE(ekb_tmp)
 ABI_DEALLOCATE(ffspl_tmp)

!DEBUG
!write(std_out,*)' psp9in : enter '
!write(std_out,*)' psp9in : indlmn(1:6,jj)'
!do jj=1,lmnmax
!write(std_out,*)indlmn(1:6,jj)
!end do
!ENDDEBUG

 ABI_DEALLOCATE( vpspll   )
 ABI_DEALLOCATE( rad      )
 ABI_DEALLOCATE( vloc     )
 ABI_DEALLOCATE( wfll     )
 ABI_DEALLOCATE( rfhi     )
 ABI_DEALLOCATE( ratm     )
 ABI_DEALLOCATE( vps      )

 ABI_DEALLOCATE( wfllatm  )
 ABI_DEALLOCATE( vlocal   )
 ABI_DEALLOCATE( chcore   )

#else
!Initialize some arguments, for portability at compile time
 indlmn=0 ; mmax=0 ; nproj=0
 ekb=zero ; epsatm=zero ; ffspl=zero ; qchrg=zero ; vlspl=zero ; xcccrc=zero ; xccc1d=zero
#endif

end subroutine psp9in
!!***
