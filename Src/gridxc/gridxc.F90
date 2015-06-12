!!@LICENSE

!******************************************************************************
! MODULE gridXC
!------------------------------------------------------------------------------
! Provides the following main XC routines:
!   atomXC   ! XC for a spherical charge distribution
!   cellXC   ! XC for a periodic unit cell
!   getXC    ! Returns the XC functional(s) being used
!   setXC    ! Sets XC functional(s) to be used by atomXC and/or cellXC
!
! Real kinds (precision) of arguments to call atomxc and cellxc
!   gridXC_std_p  ! Standard real-kind (double) precision
!   gridXC_grid_p ! Precision for grid arrays to call cellxc

! Secondary entry points for testers and lower-level programming
!   ldaxc    ! LDA-XC functionals
!   ggaxc    ! GGA-XC functionals

! Extra utilities placed here for non-siesta users
! See correspondig modules for usage documentation
!   nfft                 ! Get allowed sizes for FFTs
!   setDebugOutputUnit   ! Initialize debug report
!   closeDebugOutputFile ! Print debug report
!   myMeshBox            ! Get my processor mesh box
!   setMeshDistr         ! Set a distribution of mesh points over processors

! --------- DEPENDENCIES ------------------------------------------------------
! Modules used and what they provide:
!   alloc      : (re)allocation utility routines
!   cellsubs   : routines to find cell volume and reciprocal vectors
!   debugXC    : routines to set and print a report for debugging XC library
!   interpolation: One dimensional cubic spline interpolation routines
!   m_fft_gpfa : One dimensional Fourier transform
!   fftr       : Three dimensional Fourier transform of real functions
!   m_atomXC   : atomXC routine
!   m_bessph   : spherical Bessel function routine
!   m_cellXC   : cellXC routine
!   m_chkgmx   : routine to find the planewave cutoff of the mesh
!   m_fft3d    : Three dimensional Fourier transform of complex functions
!   m_ggaxc    : routines for GGA XC functionals
!   m_io       : routine to find and reserve an available unit for I/O
!   m_ldaxc    : routines for LDA XC functionals
!   m_minvec   : routine to find the basis of shortest lattice vectors
!   m_radfft   : radial fast Fourier transform
!   m_vdwxc    : routines for the Van der Waals functional
!   m_vv_vdwxc : routines for the Vydrov-VanVoorhis VdW functional
!   mesh1D     : utilities to manipulate 1D meshes
!   mesh3D     : routines to handle mesh arrays distributed among processors
!   moreParallelSubs : utility routines to simplify some MPI communications
!   parallel   : set and keep parameters like node index and number of nodes
!   precision  : parameters 'dp' and 'grid_p' for real kinds
!   sorting    : utility routines to sort an array of integer or real numbers
!   sys        : stopping subroutine 'die', which might be provided as an
!                external routine 
!   xcmod      : setxc routine
!
! Additional modules used in parallel compilation:
!   mpi_siesta : f90 interfaces to MPI routines, wrapped by CPU-timing calls
!
!------------- COMPILATION ----------------------------------------------------
!   A Makefile is provided in the gridXC directory to make the
! gridXC library. 
!   For a parallel compilation, compile with -DMPI
!   For CPU-time profiling of MPI calls, compile with -DMPI_TIMING
!   Parameters dp and grid_p in module precision determine the real
! kinds used in the input-output and internal variables and arrays. 
!   By default, grid_p, the kind of arrays Dens and Vxc in cellXC, is
! set to single precision, but it can be turned to double precision
! by compiling with -DGRID_DP
!   By default, dp is set to double precision. To change it, edit precision.F
!
!******************************************************************************
!
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
! Energy unit depending of internal parameter Eunit
! ------------------------ USAGE ----------------------------------------------
! You must call setXC before calling atomXC for the first time.
! A typical call sequence is:
!
!   use precision, only: dp
!   use gridXC,  only: setXC, atomXC
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
!
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
! Energy unit depending of internal parameter EUnit
! Stress in EUnit/Bohr**3
! ------------------------ USAGE ----------------------------------------------
! You must call setXC before calling cellXC for the first time.
! A typical serial program call is:
!
!   use gridXC, only: dp, grid_p
!   use gridXC, only: setXC, cellXC, nfft
!   integer  :: i, nMesh(3), nSpin
!   real(dp) :: cell(3,3), Dc, Dx, Ec, Ex, stress(3,3), 
!   real(grid_p),allocatable :: dens(:,:,:,:), Vxc(:,:,:,:)
!     Find nSpin, cell(:,:), and nMesh(:)
!   do i = 1,3
!     call nfft( nMesh(i) )   ! Increase nMesh if necessary for FFTs
!   end do
!   allocate( dens(nMesh(1),nMesh(2),nMesh(3),nSpin), &
!              Vxc(nMesh(1),nMesh(2),nMesh(3),nSpin)) )
!     Find dens(:,:,:,:) at all mesh points
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call cellXC( 0, cell, nMesh, 1,nMesh(1), 1,nMesh(2), 1,nMesh(3), &
!                nSpin, dens, Ex, Ec, Dx, Dc, stress, Vxc )
!
! A typical parallel program call is:
!
!   use gridXC, only: dp, grid_p
!   use gridXC, only: setXC, cellXC, nfft
!   integer  :: i, iSpin, myBox(2,3), nMesh(3), nSpin
!   real(dp) :: cell(3,3), Dc, Dx, Ec, Ex, stress(3,3), 
!   real(grid_p),allocatable :: dens(:,:,:,:), Vxc(:,:,:,:)
!     Find nSpin, cell(:,:), nMesh(:)
!   do i = 1,3
!     call nfft( nMesh(i) )   ! Increase nMesh if necessary for FFTs
!   end do
!     Find myBox(:,:), i.e. distribute mesh points over processors
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
!                nSpin, dens, Ex, Ec, Dx, Dc, stress, Vxc )
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
!
! subroutine setXC( n, func, auth, wx, wc )
! -----------------------------------------------------------------------------
! Sets the xc functional(s) to be used by atomXC and/or cellXC
! ------------------------- INPUT ---------------------------------------------
!     integer,         :: n       ! Number of functionals
!     character(len=*) :: func(n) ! Functional name labels
!     character(len=*) :: auth(n) ! Functional author labels
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
!           'AM05' => GGA Mattsson & Armiento, PRB, 79, 155101 (2009)
!      'PBEJsJrLO' => GGA Reparametrizations of the PBE functional by
!     'PBEJsJrHEG' => GGA   L.S.Pedroza et al, PRB 79, 201106 (2009) and
!      'PBEGcGxLO' => GGA   M.M.Odashima et al, JCTC 5, 798 (2009)
!     'PBEGcGxHEG' => GGA using 4 different combinations of criteria
!          'DRSLL' => VDW Dion et al, PRL 92, 246401 (2004)
!          'LMKLL' => VDW K.Lee et al, PRB 82, 081101 (2010)
!            'KBM' => VDW optB88-vdW of J.Klimes et al, JPCM 22, 022201 (2010)
!            'C09' => VDW V.R. Cooper, PRB 81, 161104 (2010)
!             'BH' => VDW K. Berland and Per Hyldgaard, PRB 89, 035412 (2014)
!             'VV' => VDW Vydrov-VanVoorhis, JCP 133, 244103 (2010)
!
! ------------------------ USAGE ----------------------------------------------
!   use gridXC, only: setXC
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )

! --------------------- BEHAVIOUR ---------------------------------------------
! - Stops with an error message if n is larger than internal parameter maxFunc
!
!******************************************************************************
!
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
!   use gridXC,  only: getXC
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
!
! SUBROUTINE LDAXC( AUTHOR, IREL, nSpin, D, 
!                   EPSX, EPSC, VX, VC, DVXDN, DVCDN)
! -----------------------------------------------------------------------------
! Finds the local exchange and correlation energies per electron, and their
! derivatives with respect to density and density gradient, in the
! Generalized Gradient Approximation.
! ------------------------- INPUT ---------------------------------------------
!   character(len=*) AUTHOR  ! GGA flavour (author initials)
!   integer   IREL           ! Relativistic exchange? 0=no, 1=yes
!   integer   nSpin          ! Number of spin components
!   real(dp)  D(nSpin)       ! Local electron (spin) density
!
! Allowed author values:
!     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
!           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992)
!
! ------------------------- OUTPUT --------------------------------------------
!   real(dp) EPSX            ! Exchange energy per electron
!   real(dp) EPSC            ! Correlation energy per electron
!   real(dp) VX(nspin)       ! Exchange potential, defined as dEx/dD(ispin)
!   real(dp) VC(nspin)       ! Correlation potential, defined as dEc/dD(ispin)
!   real(dp) DVXDN(nspin,nspin) ! Derivative of exchange potential with
!                                 respect the charge density, defined 
!                                 as DVx(spin1)/Dn(spin2)
!   real(dp) DVCDN(nspin,nspin) ! Derivative of correlation potential
!                                 respect the charge density, defined 
!                                 as DVc(spin1)/Dn(spin2)
! ------------------------ UNITS ----------------------------------------------
! Distances in atomic units (Bohr).
! Densities in atomic units (electrons/Bohr**3)
! Energies in Hartrees
!
!******************************************************************************
!
! SUBROUTINE GGAXC( AUTHOR, IREL, nSpin, D, GD,
!                   EPSX, EPSC, dEXdD, dECdD, dEXdGD, dECdGD )
! -----------------------------------------------------------------------------
! Finds the local exchange and correlation energies per electron, and their
! derivatives with respect to density and density gradient, in the
! Generalized Gradient Approximation.
! ------------------------- INPUT ---------------------------------------------
!   character(len=*) AUTHOR  ! GGA flavour (author initials)
!   integer   IREL           ! Relativistic exchange? 0=no, 1=yes
!   integer   nSpin          ! Number of spin components
!   real(dp)  D(nSpin)       ! Local electron (spin) density
!   real(dp)  GD(3,nSpin)    ! Gradient of electron density
!
! Allowed author values:
!           'PW91' => GGA Perdew & Wang, JCP, 100, 1290 (1994) 
!            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
!           'RPBE' => GGA Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
!         'revPBE' => GGA Zhang & Yang, PRL 80,890(1998)
!            'LYP' => GGA Becke-Lee-Yang-Parr (see subroutine blypxc)
!             'WC' => GGA Wu-Cohen (see subroutine wcxc)
!         'PBESOL' => GGA Perdew et al, PRL, 100, 136406 (2008)
!           'AM05' => GGA Mattsson & Armiento, PRB, 79, 155101 (2009)
!
! ------------------------- OUTPUT --------------------------------------------
!   real(dp) EPSX            ! Exchange energy per electron
!   real(dp) EPSC            ! Correlation energy per electron
!   real(dp) dEXdD(nSpin)    ! dEx/dDens, Ex=Int(dens*epsX)
!   real(dp) dECdD(nSpin)    ! dEc/dDens
!   real(dp) dEXdGD(3,nSpin) ! dEx/dGrad(Dens)
!   real(dp) dECdGD(3,nSpin) ! dEc/dGrad(Dens)
! ------------------------ UNITS ----------------------------------------------
! Distances in atomic units (Bohr).
! Densities in atomic units (electrons/Bohr**3)
! Energies in Hartrees
!
!******************************************************************************
!
! subroutine nfft( n )
! -----------------------------------------------------------------------------
! Changes n to the next integer allowed by the FFT routines used in cellXC
! --------------------- OPTIONAL OUTPUT ---------------------------------------
! integer:: n    ! Candidate FFT vector size, increased to
!                ! the next allowed value, if necessary
! ------------------------ USAGE ----------------------------------------------
! See usage of cellXC
!
!******************************************************************************
!
! subroutine setDebugOutputUnit()
! --------------------------------------------------------------------
! Sets debug output unit and opens file for debug output
! ------------------------ USAGE -------------------------------------
! The following call should be previous to that of atomxc and cellxc:
!   use gridXC, only: setDebugOutputUnit
!   call setDebugOutputUnit()
!
!******************************************************************************
!
! subroutine closeDebugOutputFile()
! --------------------------------------------------------------------
! Closes debug output file and copies it to node 0
! The file name is debugXC.out in serial execution and debugXC.node$$
! in parallel execution.
! ------------------------ USAGE -------------------------------------
! The following call should be made after last call to atomxc and cellxc:
!   use gridXC, only: closeDebugOutputUnit
!   call closeDebugOutputUnit()
!******************************************************************************
!
! SUBROUTINE setMeshDistr( distrID, nMesh, box, firstNode, nNodes, &
!                          nNodesX, nNodesY, nNodesZ, nBlock, &
!                          wlDistr, workload )
!
! Defines a parallel distribution of mesh points over processor nodes.
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3)  : Mesh divisions in each axis (including siesta "subpoints")
!------------------------- OPTIONAL INPUT -------------------------------------
! integer box(2,3)  : Mesh box of my processor node: 
!                     box(1,iAxis)=lower box limits, in range (0:nMesh(iAxis)-1)
!                     box(2,iAxis)=upper box limits, in same range
! integer firstNode : First node in the mesh distr.
! integer nNodes    : Total nodes in the mesh distr
! integer nNodesX   : Nodes in the X (first) axis.
!                     Must be present if present(nNodesYZ)
! integer nNodesY   : Nodes in the Y (second) axis.
!                     Must be present if present(nNodesZ)
! integer nNodesZ   : Nodes in the Z (third) axis
! integer nBlock    : Size of blocks of mesh points, in each axis, which are
!                     not splitted in different nodes. It must be a factor of
!                     all of nMesh(1:3). If box is also present, nBlock must be
!                     a factor of all box(1,1:3) and box(2,1:3)+1.
!                     nBlock corresponds to nsm (lateral size of "superpoints")
!                     in siesta mesh terminology.
! integer wlDistr   : Distr. index of workload array
! real(gp) workload(0:,0:,0:) ! Approx. relative workload of mesh points. 
!                     Must be nonnegative at all points and have nonzero sum.
!-------------------------- INPUT and OUTPUT ----------------------------------
! integer distrID : ID assigned to the mesh distrib.
!----------------------------- USAGE ------------------------------------------
!    Arguments box, firstNode, nNodes, nNodesXYZ, and nBlock are provided to
! force compatibility with user distributions, but they should be avoided
! otherwise, since they difficult the optimal distribution of mesh points.
!
!    Typical usage to create a uniform 3D mesh distribution:
! use gridXC, only: setMeshDistr
! integer:: nMesh(3)
! integer,save:: myDistr=-1
! ... Find nMesh
! call setMeshDistr( myDistr, nMesh )
!
!    Typical usage to distribute evenly the workload of each node:
! use gridXC, only: setMeshDistr
! integer:: nMesh(3)
! integer,save:: newDistr=-1, oldDistr=-1
! real(gp),allocatable:: workload(:,:,:)
! ... Find nMesh
! call setMeshDistr( oldDistr, nMesh )
! call myMeshBox( oldDistr, nMesh, box )
! allocate( workload(box(1,1):box(2,1), &
!                    box(1,2):box(2,2), &
!                    box(1,3):box(2,3)) )
! ... Find approximate CPU workload associated to each point of box
! call setMeshDistr( newDistr, nMesh, wlDistr=oldDistr, workload=workload )
!---------------------------- BEHAVIOUR ---------------------------------------
! In serial execution (totNodes==1) it simply returns distrID=0, irrespective
!   of all arguments. Parameter totNodes is obtained from module parallel.
! If the input distribution ID is still valid (i.e. consistent with the other
!   input arguments), the same value is returned. If it points to an existing
!   distribution that is no longer consistent, the old distribution ID is 
!   freed before returning with a new distrID. This makes it convenient to 
!   make succesive calls with the same distrID but different other arguments
!   (e.g. different workloads in different iterations).
! New IDs are never repeated, even if they identify the same distribution.
! If box is present, all other optional arguments are ignored. The different
!   node mesh boxes should be a nonoverlapping partition of all the mesh 
!   points, but this is NOT checked.
! If firstNode is not present, firstNode=0 is assumed.
! If nNodes is not present all nodes are used (nNodes=totNodes).
! It stops with an error message in the following cases:
!   - nNodes is present and it is smaller than 1 or larger than totNodes.
!   - nNodesZ is present but either nNodesX or nNodesY are not present.
!   - nNodesY is present but nNodesX is not present.
!   - nBlock is present and it is not a factor of one of nMesh(1:3).
!   - box and nBlock are present, and nBlock is not a factor of any of 
!     box(1,1:3) or box(2,1:3)+1.
!   - workload is present but wlDistr is not present.
!   - workload is present and negative at any mesh point, or zero at all 
!     points of the unit cell mesh.
! If nNodesXYZ are not present, the distribution of nodes over each axis is
!   optimized, in the sense of leading to node mesh boxes as cubic as possible. 
!   If only nNodesX is present, it is optimized over the Y and Z axes.
! If nBlock is not present, nBlock=1 is assumed.
! If workload is present, its distribution over processors is optimized, in the
!   sense of load balance. If it is not present, the even distribution of mesh 
!   points is optimized (equivalent to using the same workload for all points).
!--------------------------- ALGORITHMS ---------------------------------------
! For uniform distributions (workload not present):
! - nNodes is first fatorized in its prime factors.
! - All possible distributions of (products of) factors over axes are tried,
!   and that with lowest dispersion of nMesh(axis)/factor(axis) is selected.
!   This is done only over axes not constrained by nNodesXYZ, if present, and
!   respecting that axis factors must be multiples of nBlock, if present.
! - nMesh(axis) is divided in nNodes(axis) boxes. If there is a rest, the first
!   rest(axis) of the nNodes(axis) are given an extra point. Again, this is 
!   done with blocks of nBlock points (if present) rather than single points.
!
! For nonuniform distributions (workload present):
! - nNodes is first fatorized in its prime factors.
! - Beginning with a box of the entire cell size, each box is divided in 
!   nParts=factor (in order of decreasing factors). The division axis is that
!   along which the spatial dispersion of the box workload is maximum. The 
!   division surfaces are chosen to split the workload as evenly as possible.
!   If the division surface is in a void region (with zero workload), it is
!   placed at the midpoint between the nonzero workload regions. 
! - The divided boxes are assigned to a contiguous set of 
!   boxNodes=boxNodes/nParts. This process is repeated until boxNodes=1.
!
! A stored distribution is defined by nMesh and its node boxes. There may be up
!   to maxDistr distributions defined simultaneously. Each distribution may be
!   identified simultaneously by up to maxDistrID identifiers (for example by 
!   different calling routines), but its information is stored only once. 
!   This is ensured by comparing the boxes of a new distrID with those of all
!   previously defined distributions. If they coincide, the distrID is simply
!   added to the existing distribution.
! Each distrID is never repeated. When a distrID is erased (freed), the
!   distribution itself is not erased, unless all its IDs have been erased.
!   This prevents that a calling routine may erase a distribution that is
!   still being used by other routines.
!
!******************************************************************************
!
! SUBROUTINE myMeshBox( nMesh, distrID, box )

! Finds the mesh box of the local processor in a parallel mesh distribution.
! Equivalent to nodeMeshBox with node=myNode
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3) : Mesh divisions in each axis
! integer distrID  : Mesh-distribution ID
!----------------------------- OUTPUT -----------------------------------------
! integer box(2,3) : Mesh box: box(1,:)=lower bounds, box(2,:)=upper bounds
!----------------------------- USAGE ------------------------------------------
!   Typical usage to find my box of mesh points:
! use gridXC, only: myMeshBox
! integer:: myBox(2,3), nMesh(3)
! integer,save:: myDistr=-1
! ...Find nMesh and myDistr
! call myMeshBox( nMesh, myDistr, myBox )
!---------------------------- BEHAVIOUR ---------------------------------------
! If iDistr=0, it returns box(1,:)=0 and  box(2,:)=nMesh(:)-1 
! If node stores no mesh points, its (empty) box returns with box(1,:)>box(2,:)
!
!******************************************************************************


MODULE gridXC

! Real kinds (precision) of arguments

  USE precision, only: dp      ! Standard real-kind (double) precision
  USE precision, only: grid_p  ! Precision for grid arrays

! Main entry routines of gridXC library
!--------------------------------------------------------------
  use gridxc_config, only: gridxc_init  ! initialization

  ! XC for a spherical charge distribution
  USE m_atomXC, only: gridxc_atomXC => atomXC 

  USE m_cellXC, only: gridxc_cellXC => cellXC  ! XC for a periodic unit cell
  USE xcmod,    only: gridxc_getXC => getXC   ! Returns XC functional(s)
  USE xcmod,    only: gridxc_setXC => setXC   ! Sets XC functional(s)

! Secondary entry points for testers and lower-level programming
  USE m_ldaxc,  only: gridxc_ldaxc => ldaxc    ! LDA-XC functionals
  USE m_ggaxc,  only: gridxc_ggaxc => ggaxc    ! GGA-XC functionals

! Extra utilities placed here for non-siesta users
! See correspondig modules for usage documentation
  USE m_fft_gpfa,    only: gridxc_nfft => nfft  ! Get allowed sizes for FFTs
!----------------------------------------------------------------------
!-----------------------------------
#ifdef DEBUG_XC
  USE debugXC,  only: setDebugOutputUnit   ! Set debug report
  USE debugXC,  only: closeDebugOutputFile ! Print debug report
#endif
  USE mesh3d,   only: myMeshBox            ! Get my processor mesh box
  USE mesh3d,   only: setMeshDistr         ! Set a distribution of mesh
                                           ! points over parallel processors
  PUBLIC

END MODULE gridXC
  
