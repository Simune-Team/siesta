!!@LICENSE

!******************************************************************************
! MODULE siestaXC
!------------------------------------------------------------------------------
! Provides the following routines:
!   atomXC   ! XC for a spherical charge distribution
!   cellXC   ! XC for a periodic unit cell
!   getXC    ! Returns the XC functional(s) being used
!   setXC    ! Sets XC functional(s) to be used by atomXC and/or cellXC
!
! --------- DEPENDENCIES ------------------------------------------------------
! Routines called in serial: 
!   meshKcut, RECLAT, TIMER, VOLCEL
! Modules used:
!   alloc     : provides (re)allocation utility routines
!   m_atomXC  : provides atomXC routine
!   m_cellXC  : provides cellXC routine
!   m_ggaxc   : provides routines for GGA XC functionals
!   m_ldaxc   : provides routines for LDA XC functionals
!   m_radfft  : provides radial fast Fourier transform
!   m_vdwxc   : povides routines for the Van der Waals functional
!   mesh1D    : provides utilities to manipulate 1D meshes
!   mesh3D    : provides routines to handle mesh arrays distributed
!               among processors
!   precision : defines parameters 'dp' and 'grid_p' for real kinds
!   sys       : provides the stopping subroutine 'die'
!   xcmod     : provides setxc routine
!   TO BE COMPLETED...
! Additional modules used in parallel:
!   mpi_siesta
!
!------------- COMPILATION ----------------------------------------------------
!   A Makefile is provided in the siestaXC directory to make the
! siestaXC library. 
!   For a parallel compilation, compile with -DMPI
!   Parameters dp and grid_p in module precision determine the real
! kinds used in the input-output and internal variables and arrays. 
! By default, grid_p, the kind of arrays Dens and Vxc in cellXC, is
! set to single precision, but it can be turned to double precision
! by compiling with -DGRID_DP
!
!------------- USAGE ----------------------------------------------------------
! You must call setXC before calling cellXC for the first time.
! See usage examples in the headers of each routine below.
!
!******************************************************************************
! subroutine atomXC( irel, nr, maxr, rmesh, nSpin, Dens, Ex, Ec, Dx, Dc, Vxc )
!------------------------------------------------------------------------------
! Finds total exchange-correlation energy and potential for a
! spherical electron density distribution.
! This version implements the Local (spin) Density Approximation and
! the Generalized-Gradient-Aproximation with the 'explicit mesh 
! functional' approach of White & Bird, PRB 50, 4954 (1994).
! Gradients are 'defined' by numerical derivatives, using 2*nn+1 mesh
!   points, where nn is a parameter defined below
! Ref: L.C.Balbas et al, PRB 64, 165110 (2001)
! Written by J.M.Soler using algorithms developed by 
!   L.C.Balbas, J.L.Martins and J.M.Soler, Dec.1996
! Van der Waals functional added by J.M.Soler, Jul.2008, as explained in
!   G.Roman-Perez and J.M.Soler, PRL 103, 096102 (2009)
! ------------------------- INPUT ---------------------------------------------
! INTEGER  irel         : Relativistic exchange? (0=>no, 1=>yes)
! INTEGER  nr           : Number of radial mesh points
! INTEGER  maxr         : Physical first dimension of Dens and Vxc
! REAL(dp) rmesh(nr)    : Radial mesh points. Must be nr.le.maxr
! INTEGER  nSpin        : nSpin=1 => unpolarized; nSpin=2 => polarized
! REAL(dp) Dens(maxr,nSpin) : Total (nSpin=1) or spin (nSpin=2) electron
!                            density at mesh points
! ------------------------- OUTPUT --------------------------------------------
! REAL(dp) Ex              : Total exchange energy
! REAL(dp) Ec              : Total correlation energy
! REAL(dp) Dx              : IntegralOf( rho * (eps_x - v_x) )
! REAL(dp) Dc              : IntegralOf( rho * (eps_c - v_c) )
! REAL(dp) Vxc(maxr,nSpin) : (Spin) exch-corr potential
! ------------------------ UNITS ----------------------------------------------
! Distances in atomic units (Bohr).
! Densities in atomic units (electrons/Bohr**3)
! Energy unit depending of parameter Eunit below
! ------------------------ USAGE ----------------------------------------------
! You must call setXC before calling atomXC for the first time.
! A typical call sequence is:
!
!   use precision, only: dp
!   use siestaXC,  only: setXC, atomXC
!   integer  :: nr, nSpin
!   real(dp) :: Dc, Dx, Ec, Ex
!   real(dp),allocatable :: dens(:,:), rMesh(:,:), Vxc(:,:)
!     Find nr and nSpin
!   allocate( dens(nr,nSpin), rMesh(nr), Vxc(nr,nSpin) )
!     Find rMesh(:) and dens(:,:) at all mesh points
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call atomXC( 0, nr, nr, rmesh, nSpin, Dens, Ex, Ec, Dx, Dc, Vxc )
!
! --------- BEHAVIOUR ---------------------------------------------------------
! Stops and prints an error message if maxr<nr
! Stops and prints an error message if functl is not one of LDA, GGA, or VDW
!
!******************************************************************************
! subroutine cellXC( irel, cell, nMesh, lb1, ub1, lb2, ub2, lb3, ub3, 
!    .               nSpin, dens, Ex, Ec, Dx, Dc, stress, Vxc, dVxcdD )
!------------------------------------------------------------------------------
! Finds total exchange-correlation energy and potential in a
!   periodic cell.
! This version implements the Local (spin) Density Approximation and
!   the Generalized-Gradient-Aproximation with the 'explicit mesh 
!   functional' approach of White & Bird, PRB 50, 4954 (1994).
! Gradients are 'defined' by numerical derivatives, using 2*nn+1 mesh
!   points, where nn is a parameter defined below
! Ref: L.C.Balbas et al, PRB 64, 165110 (2001)
! Wrtten by J.M.Soler using algorithms developed by 
!   L.C.Balbas, J.L.Martins and J.M.Soler, Dec.1996 - Aug.1997
! Parallel version written by J.Gale. June 1999.
! Argument dVxcdD added by J.Junquera. November 2000.
! Adapted for multiple functionals in the same run by J.Gale 2005
! Van der Waals functional added by J.M.Soler, Jan.2008, as explained in
!   G.Roman-Perez and J.M.Soler, PRL 103, 096102 (2009)
! ------------------------- INPUT ---------------------------------------------
! integer  irel        : Relativistic exchange? (0=>no, 1=>yes)
! real(dp) cell(3,3)   : Unit cell vectors cell(ixyz,ivector)
! integer  nMesh(3)    : Total mesh divisions of each cell vector
! integer  lb1,lb2,lb3 : Lower bounds of arrays dens, Vxc, dVxcdD
! integer  ub1,ub2,ub3 : Upper bounds of arrays dens, Vxc, dVxcdD
! integer  nSpin       : nSpin=1 => unpolarized; nSpin=2 => polarized;
!                        nSpin>2 => non-collinear polarization
! real(grid_p) dens(lb1:ub1,lb2:ub2,lb3:ub3,nSpin) : Total (nSpin=1) or 
!                        spin (nSpin=2) electron density at mesh points
!                        For non-collinear polarization, the density
!                        matrix is given by: dens(1)=D11, dens(2)=D22,
!                        dens(3)=Real(D12), dens(4)=Im(D12)
! ------------------------- OUTPUT --------------------------------------------
! real(dp) Ex             : Total exchange energy per unit cell
! real(dp) Ec             : Total correlation energy per unit cell
! real(dp) Dx             : IntegralOf( rho * (eps_x - v_x) ) in unit cell
! real(dp) Dc             : IntegralOf( rho * (eps_c - v_c) ) in unit cell
! real(dp) stress(3,3)    : xc contribution to the stress tensor, in unit
!                           cell, assuming constant density (not charge),
!                           i.e. r->r' => rho'(r') = rho(r)
!                           For plane-wave and grid (finite diff) basis
!                           sets, density rescaling gives an extra term
!                           (not included) (Dx+Dc-Ex-Ec)/cell_volume for
!                           the diagonal elements of stress. For other
!                           basis sets, the extra term is, in general:
!                           IntegralOf(v_xc * d_rho/d_strain) / cell_vol
! real(grid_p) Vxc(lb1:ub1,lb2:ub2,lb3:ub3,nSpin) : (Spin) xc potential
! ------------------------ OPTIONAL OUTPUT ------------------------------------
! real(grid_p) dVxcdD(lb1:ub1,lb2:ub2,lb3:ub3,nSpin*nSpin) : Derivatives
!                           of xc potential respect to charge density
!                           Available only for LDA
! ------------------------ UNITS ----------------------------------------------
! Distances in atomic units (Bohr).
! Densities in atomic units (electrons/Bohr**3)
! Energy unit depending of parameter EUnit below
! Stress in EUnit/Bohr**3
! ------------------------ USAGE ----------------------------------------------
! You must call setXC before calling cellXC for the first time.
! A typical serial program call is:
!
!   use precision, only: dp, grid_p
!   use siestaXC,  only: setXC, cellXC
!   integer  :: nMesh(3), nSpin
!   real(dp) :: cell(3,3), Dc, Dx, Ec, Ex, stress(3,3), 
!   real(grid_p),allocatable :: dens(:,:,:,:), Vxc(:,:,:,:)
!     Find nSpin, cell(:,:), and nMesh(:)
!   allocate( dens(nMesh(1),nMesh(2),nMesh(3),nSpin), &
!              Vxc(nMesh(1),nMesh(2),nMesh(3),nSpin)) )
!     Find dens(:,:,:,:) at all mesh points
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call cellXC( 0, cell, nMesh, 1,nMesh(1), 1,nMesh(2), 1,nMesh(3), &
!                nSpin, dens, Ex, Ex, Dx, Dc, stress, Vxc )
!
! A typical parallel program call is:
!
!   use precision, only: dp, grid_p
!   use siestaXC,  only: setXC, cellXC
!   integer  :: iSpin, myBox(2,3), nMesh(3), nSpin
!   real(dp) :: cell(3,3), Dc, Dx, Ec, Ex, stress(3,3), 
!   real(grid_p),allocatable :: dens(:,:,:,:), Vxc(:,:,:,:)
!     Find nSpin, cell(:,:), nMesh(:), and myBox(:,:)
!   allocate( dens(myBox(1,1):myBox(2,1),        &
!                  myBox(1,2):myBox(2,2),        &
!                  myBox(1,3):myBox(2,3),nSpin), &
!              Vxc(myBox(1,1):myBox(2,1),        &
!                  myBox(1,2):myBox(2,2),        &
!                  myBox(1,3):myBox(2,3),nSpin) )
!   do i3 = myBox(1,3),myBox(2,3)
!   do i2 = myBox(1,2),myBox(2,2)
!   do i1 = myBox(1,1),myBox(2,1)
!     do iSpin = 1,nSpin
!       dens(i1,i2,i3,iSpin) = (spin)density at point (i1,i2,i3)
!     end do
!   end do
!   end do
!   end do
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call cellXC( 0, cell, nMesh, myBox(1,1), myBox(2,1), &
!                                myBox(1,2), myBox(2,2), &
!                                myBox(1,3), myBox(2,3), &
!                nSpin, dens, Ex, Ex, Dx, Dc, stress, Vxc )
!
! IMPORTANT: arrays dens, Vxc, and dVxcdD may be alternatively 
! allocated and initialized with indexes beginning in 0 or 1, 
! or even use a single-index array, e.g.:
!   real(grid_p),allocatable :: dens(:,:), Vxc(:,:)
!   myMesh(1:3) = myBox(2,:) - myBox(1,:) + 1
!   myPoints = myMesh(1)*myMesh(2)*myMesh(3)
!   allocate( dens(myPoints,nSpin), Vxc(myPoints,nSpin) )
!   iPoint = 0
!   do i3 = myBox(1,3),myBox(2,3)
!   do i2 = myBox(1,2),myBox(2,2)
!   do i1 = myBox(1,1),myBox(2,1)
!     iPoint = iPoint + 1
!     do iSpin = 1,nSpin
!       dens(iPoint,iSpin) = (spin)density at point (i1,i2,i3)
!     end do
!   end do
!   end do
!   end do
! but the call to cellXC must still be as given above, i.e. the
! arguments lb1,ub1,... must be the lower and upper bounds of the
! mesh box stored by each processor (not the allocated array bounds).
! However, the arrays size MUST be (ub1-lb1+1)*(ub2-lb2+1)*(ub3-lb3+1)
! The processor mesh boxes must not overlap, and they must cover the 
! entire unit cell mesh. This is NOT checked inside cellXC.
!
! --------- BEHAVIOUR ---------------------------------------------------------
! - Stops and prints a warning if functl is not one of LDA, GGA, or VDW
! - The output values of Ex, Ec, Dx, Dc, and stress, are integrals over
!   the whole unit cell, not over the mesh box of the local processor
! - Since the exchange and correlation part is usually a small fraction
!   of a typical electronic structure calculation, this routine has
!   been coded with emphasis on simplicity and functionality, not in
!   efficiency.
!
!******************************************************************************
! subroutine setXC( n, func, auth, wx, wc )
! -----------------------------------------------------------------------------
! Sets the xc functional(s) to be used by atomXC and/or cellXC
! ------------------------- INPUT ---------------------------------------------
!     integer,         :: n       ! Number of functionals
!     character(len=*),:: func(n) ! Functional name labels
!     character(len=*),:: auth(n) ! Functional author labels
!     real(dp),        :: wx(n)   ! Functional weights for exchange
!     real(dp),        :: wc(n)   ! Functional weights for correlation
!
! Allowed functional/author values:
! XCfunc: 
!   'LDA' or 'LSD' => Local density approximation
!            'GGA' => Generalized gradients approx.
!            'VDW' => Van der Waals functional
! XCauth:
!     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
!           'PW91' => GGA Perdew & Wang, JCP, 100, 1290 (1994) 
!           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
!                     the local density limit of the next:
!            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
!           'RPBE' => GGA Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
!         'revPBE' => GGA Zhang & Yang, PRL 80,890(1998)
!            'LYP' => GGA Becke-Lee-Yang-Parr (see subroutine blypxc)
!             'WC' => GGA Wu-Cohen (see subroutine wcxc)
!         'PBESOL' => GGA Perdew et al, PRL, 100, 136406 (2008)
!          'DRSLL' => VDW Dion et al, PRL 92, 246401 (2004)

! ------------------------ USAGE ----------------------------------------------
!   use siestaXC, only: setXC
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )

! --------------------- BEHAVIOUR ---------------------------------------------
! - Stops with an error message if n is larger than internal parameter maxFunc
!
!******************************************************************************
! subroutine getXC( n, func, auth, wx, wc )
! -----------------------------------------------------------------------------
! Returns the xc functional(s) that has been previously set
! --------------------- OPTIONAL OUTPUT ---------------------------------------
!     integer         :: n       ! Number of functionals
!     character(len=*):: func(n) ! Functional name labels
!     character(len=*):: auth(n) ! Functional author labels
!     real(dp)        :: wx(n)   ! Functional weights for exchange
!     real(dp)        :: wc(n)   ! Functional weights for correlation
!
! ------------------------ USAGE ----------------------------------------------
!   use precision, only: dp
!   use siestaXC,  only: getXC
!   integer,parameter:: maxFunc = 10
!   character(len=20):: func(maxFunc), auth(maxFunc)
!   real(dp):: wx(maxFunc), wc(maxFunc)
!   call getXC( n, func, auth, wx, wc )
!
! --------------------- BEHAVIOUR ---------------------------------------------
! - Stops with an error message if called before setXC
! - Does not change any output array whose size is smaller than nFunc
!
!******************************************************************************

MODULE siestaXC

! Real kinds (precision) of arguments
  USE precision, only: siestaXC_std_p  => dp     ! Standard real-kind precision
  USE precision, only: siestaXC_grid_p => grid_p ! Precision for grid arrays

! Main entry routines of siestaXC library
  USE m_atomXC, only: atomXC   ! XC for a spherical charge distribution
  USE m_cellXC, only: cellXC   ! XC for a periodic unit cell
  USE xcmod,    only: getXC    ! Returns XC functional(s)
  USE xcmod,    only: setXC    ! Sets XC functional(s)

! Secondary entry points for testers and lower-level programming
  USE m_ldaxc,  only: ldaxc    ! LDA-XC functionals
  USE m_ggaxc,  only: ggaxc    ! GGA-XC functionals

! Extra utilities placed here for non-siesta users
! See correspondig modules for usage documentation
  USE fft1d,    only: nfft                 ! Get allowed sizes for FFTs
  USE alloc,    only: alloc_report         ! Set and print allocation report
  USE debugXC,  only: setDebugOutputUnit   ! Set debug report
  USE debugXC,  only: closeDebugOutputFile ! Print debug report
  USE m_timer,  only: timer_report         ! Print CPU time report
  USE mesh3d,   only: myMeshBox            ! Get my processor mesh box
  USE mesh3d,   only: setMeshDistr         ! Set a distribution of mesh
                                           ! points over parallel processors
  PUBLIC

END MODULE siestaXC
  
