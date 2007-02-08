	De: 	  jose.soler@uam.es
	Asunto: 	m_siesta comentado
	Fecha: 	26 de enero de 2007 17:14:19 GMT+01:00
	Para: 	  albertog@icmab.es, josem.cela@bsc.es, manuel.quero@bsc.es

Alberto:

Aqui va el prometido modulo m_siesta con las variables
comentadas. Espero que ayude algo a Manuel a entender
el codigo, y tambien a distribuir racionalmente las
variables en modulos mas pequeños.

No he incluido el tipo integer(is). La verdad es que
a mi no me parece muy necesario, pues el integer
estandard nos ha bastado siempre en todas las maquinas.
Pero si os parece necesario, pues vosotros mismos.
Unicamente hay que tener cuidado porque el nombre 'is'
se usa como indice de especies. Podria cambiarse este
por ispec, o algo similar.

Jose
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2007.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
!----------------------------------------------------------------------------
!
! Provisional module to contain the variables and arrays of the main
! siesta program, vesion 2.1. Written by J.M.Soler. Jan.2007.
!
! Frequently used abreviations:
!  BZ : Brilouin Zone (Unit cell in reciprocal space)
!  CG : Conjugate Gradient (minimization method)
!  DM : Density Matrix
!  IO : Input/Output
!  k  : Reciprocal-space vector
!  H  : Hamiltonian matrix
!  KB : Kleinman-Bylander (nonlocal pseudopotential)
!  MD : Molecular Dynamics
!  PR : Parrinello-Rahman (variable-cell dynamics)
!  S  : Overlap matrix
!  SCF: Self Consistent Field (between electron density and potential)
!
!----------------------------------------------------------------------------

MODULE m_siesta

use files,     only: label_length ! Max. characters of file-name labels
use parse,     only: parsed_line  ! Derived type with contents of a line read
use precision, only: dp           ! Double-precision real type

implicit none

PUBLIC

integer:: broyden_maxit ! Max. iterations in Broyden geometry relaxation
integer:: fincoor ! Final step in geommetry iteration
integer:: i       ! General-purpose index
integer:: ia      ! Atom index
integer:: ia1     ! Atom index
integer:: ia2     ! Atom index
integer:: iadispl ! Displaced atom index
integer:: ianneal ! Annealing option, read in redata and passed to anneal
integer:: idyn    ! Geommetry relaxation/dynamics option
integer:: ifa     ! Calculate forces? 0=>no, 1=>yes
integer:: ifinal  ! Last geommetry iteration step, for some types of dynamics
integer:: ihmat   ! Calculate hamiltonian? 0=>no, 1=>yes
integer:: ihuge=1073741823 ! Huge integer
integer:: iiscf   ! Auxiliary to call pulayx, generally the current SCF iter
integer:: ik      ! k-vector index
integer:: ind     ! Index of sparse matrix elements
integer:: inicoor ! Initial step in geommetry iteration
integer:: io      ! Basis-orbital index
integer:: ioptlwf ! Order-N functional option, read in redata, used in ordern
integer:: iord    ! Order of density matrix extrapolation between geom. steps
integer:: iquench ! Quenching option, read in redata, used in dynamics routines
integer:: is      ! Species index
integer:: iscf    ! SCF iteration step
integer:: isel    ! Auxiliary to call neighb: avoid double count of atom pairs?
integer:: isolve  ! Option to find density matrix: 0=>diag, 1=>order-N
integer:: ispin   ! Spin index
integer:: istart  ! First geommetry iteration step for certain types of dynamics
integer:: istep   ! Geommetry iteration step
integer:: istp    ! Geommetry iteration step starting in istp=1
integer:: istpsave! Temporary to save and restore istp when reading DM from file
integer:: istr    ! Calculate stress? 0=>no, 1=>yes
integer:: iu      ! I/O file unit
integer:: iunit   ! Physical-units option: 1=>(eV,Ang), 2=>(Ry,Bohr)
integer:: iv      ! Vector index. Not used
integer:: ix      ! Cartesian coordinate index
integer:: ixdispl ! Displaced cartesian coordinate
integer:: j       ! General-purpose index
integer:: ja      ! Atom index
integer:: jamin   ! Closest atom index
integer:: jo      ! Orbital index
integer:: jx      ! Cartesian coordinate index
integer:: kscell(3,3)=0 ! Matrix of supercell defining the BZ integration grid:
                        !    kvec = 2*pi*Inverse(cell*kscell)
integer:: level   ! Option for allocation report level of detail
integer:: maxk    ! Max. number of k vectors for BZ integration
integer:: maxkb   ! Max. number of Kleinman-Bylander proyectors
integer:: maxna=200 ! Max. number of neighbor atoms
integer:: maxnh=1 ! Max. number of nonzero H matrix elements
integer:: maxsav  ! Number of previous density matrices used in Pulay mixing
integer:: maxwk   ! Number of k-points for wavefunction printout
integer:: mscell(3,3)=(/ 1,0,0, 0,1,0, 0,0,1 /) ! Matrix of auxiliary supercell 
integer:: mullipop! Option for Mulliken population level of detail
integer:: nauxpul=1 ! Max. required size of auxiliary sparse Pulay matrices
integer:: nbcell  ! Auxiliary to call shaper routine
integer:: nbk=1   ! Number of band k points
integer:: ncgmax  ! Max. number of conjugate gradient steps in order-N minim.
integer:: neigmin ! Min. number of eigenstates (per k point)
integer:: neigwanted ! Wanted number of eigenstates (per k point)
integer:: nhist   ! Number of histogram intervals in projected density of states
integer:: nkick   ! Period between "kick" steps in SCF iteration
integer:: nkpnt=1 ! Max. required number of integration k points
integer:: nkpol   ! Number of k points for the calculation of the polarization
integer:: nh      ! Max. nonzero elements in a row of the hamiltonian matrix
integer:: nmove   ! Number of geometry iterations
integer:: nnamax=1! Max. required number of neighbor atoms
integer:: nnia    ! Number of neighbors of a given atom
integer:: no_l=1  ! Number of basis orbitals stored in one node
integer:: ns=1    ! Number of chemical species
integer:: nsc(3)=1! Diagonal elements of supercell vectors
integer:: nscf    ! Number of SCF iteration steps
integer:: nscold(3)=0 ! Previous diagonal elements of supercell vectors
integer:: nspin=1 ! Number of spin components
integer:: ntcon   ! Total number of position constraints imposed
integer:: ntm(3)  ! Number of integration mesh divisions of each cell vector
integer:: nwk     ! Number of k-points for wavefunction printout
integer:: nxij=1  ! Number of neighbors of a given basis orbital

integer, allocatable :: &
  jna(:) ! Index of neighbor atoms

integer, pointer, save :: &
  listh(:),       &! List of possibly-nonzero hamiltonian matrix elements 
  listhold(:),    &! Same as listh in previous geometry
  listhptr(:),    &! First nonzero matrix element of each matrix row
  listhptrold(:), &! Same as listhptr in previous geometry
  numh(:),        &! Number of nonzero matrix elements of each row
  numhold(:)       ! Same as numh in previous geometry

real(dp):: bcell(3,3) ! Reciprocal unit vectors
real(dp):: beta     ! Inverse temperature for Chebishev expansion.
                    ! Read in redata and used in ordern
real(dp):: bulkm    ! Bulk modulus
real(dp):: charnet  ! Net electric charge
real(dp):: cfmax    ! Max. constrained force
real(dp):: const    ! Temporary factor (used in calculation of Harris forces)
real(dp):: cstress(3,3) ! Constrained stress tensor
real(dp):: dDmax    ! Max. change in density matrix elements
real(dp):: dDtol    ! Tolerance in change of DM elements to finish SCF iteration
real(dp):: DEharr   ! Harris functional term of the total energy
real(dp):: DEna     ! Neutral-atom energy term, calculated in dnaefs
real(dp):: dipol(3) ! Electric dipole vector
real(dp):: dt       ! Time step in dynamics
real(dp):: DUext    ! Interaction energy with external electric field,
                    ! calculated in dhscf
real(dp):: DUscf    ! Electrostatic energy of (rhoscf-rhoatm), calc. in dhscf
real(dp):: Dxc      ! Integral((epsxc-Vxc)*Rho), calculated in dhscf (not used)
real(dp):: dx       ! Atomic displacement used to calculate Hessian matrix
real(dp):: dxmax    ! Max. atomic displacement allowed during geom. relaxation
real(dp):: e1=1     ! Lower bound of calculated eigenvalues
real(dp):: e2=-1    ! Upper bound of calculated eigenvalues
real(dp):: sigma    ! Energy width used to convolute partial density of states
real(dp):: Ecorrec  ! Energy term eta*DeltaQ, calculated in ordern
real(dp):: ef       ! Fermi energy
real(dp):: Eharrs   ! Harris-functional total energy
real(dp):: Eharrs1  ! Same as Eharrs, but preserved in grid cell sampling
real(dp):: Eions    ! Self-energy of isolated ions
real(dp):: Ekin     ! Kinetic energy of electrons, calculated in kinefsm
real(dp):: Ekinion  ! Kinetic energy of ions
real(dp):: Elast    ! Total energy in the previous SCF iteration
real(dp):: Emad     ! Madelung energy term, calculated in madelung
real(dp):: Emeta    ! Metadynamics energy contribution calculated in meta
real(dp):: Ena      ! Neutral-atom term in the total energy, calculated in naefs
real(dp):: Enaatm   ! Integral of Vna * rhoatm, calculated in dhscf
real(dp):: Enascf   ! Integral of Vna * rhoscf, calculated in dhscf
real(dp):: Enl      ! Non-local pseudopot. energy term, calculated in nlefsm
real(dp):: Entrop   ! Temporary to call diagon
real(dp):: Entropy  ! Entropy due to electron state occupations, calc. in diagon
real(dp):: eta(2)   ! Chemical-potential param. Read by redata, used in ordern
real(dp):: etol     ! Relative tol. in CG minim, read by redata, used in ordern
real(dp):: Etot     ! Total electronic energy
real(dp):: Etot_output ! Temporary to select among different total-energy opts
real(dp):: Exc      ! Exchange-correlation energy, calculated in dhscf
real(dp):: E0       ! Non-SCF part of total energy
real(dp):: Emm      ! Classical two-body term, calculated in twobody
real(dp):: factor   ! Temporary variable
real(dp):: fmax     ! Maximum atomic force
real(dp):: fmean    ! Mean atomic force (total force / number of atoms)
real(dp):: FreeE    ! Free energy
real(dp):: fres     ! Root mean squared atomic force
real(dp):: ftem     ! Temporary force value
real(dp):: ftol     ! Force tolerance to stop geometry relaxation
real(dp):: ftot(3)  ! Total force
real(dp):: g2cut    ! Required planewave cutoff of real-space integration mesh
real(dp):: g2max    ! Actual planewave cutoff of real-space integration mesh
real(dp):: getot    ! Total energy including thermostat dynamic variables
real(dp):: kcutof=0 ! Cutoff of the BZ integration mesh
real(dp):: kdispl(3)=0 ! Displaced origin of the BZ integration mesh
real(dp):: kn       ! Kinetic energy of the Nose' thermostat
real(dp):: kpr      ! Kinetic energy of the Parrinello-Rahman variables
real(dp):: mn       ! Mass of Nose thermostat
real(dp):: mpr      ! Mass of Parrinello-Rahman variables
real(dp):: occtol   ! Occupancy threshold to build DM
real(dp):: Pint     ! Instantaneous pressure
real(dp):: Pmol     ! "Molecular pressure" (discounting Virial term)
real(dp):: Psol     ! Solid (periodic system) pressure
real(dp):: qaux     ! Auxiliary to call spnvec (not used later)
real(dp):: qspin(4) ! Total spin populations
real(dp):: qsol     ! Total unnormalized electron charge
real(dp):: rcoor    ! Cutoff radius of Localized Wave Functions in ordern
real(dp):: rcoorcp  ! Cutoff radius to find Fermi level by projection in ordern
real(dp):: rijmin   ! Min. permited interatomic distance without warning
real(dp):: rmax_bonds ! Cutoff length for bond definition
real(dp):: rmax     ! Max. required interatomic distance
real(dp):: rmaxh    ! Max. distance between non-zero hamiltonian matrix orbitals
real(dp):: rmin     ! Min. interatomic distance
real(dp):: r2min    ! Min. interatomic distance squared
real(dp):: scell(3,3) ! Supercell unit vectors
real(dp):: stot     ! Total spin polarization
real(dp):: stress(3,3) ! Stress tensor = d_E/d_strain
real(dp):: strtol   ! Stress tolerance in relaxing the unit cell
real(dp):: svec(3)  ! Total spin vector in non-collinear-spin calculations
real(dp):: taurelax ! Relaxation time to reach desired T and P in anneal
real(dp):: tempinit ! Initial ionic temperature read in redata
real(dp):: tempion=0 ! Ionic temperature
real(dp):: tiny=1.e-15_dp ! Small tolerance used to compare k points
real(dp):: tp       ! Target pressure. Read in redata. Used in dynamics routines
real(dp):: ts       ! Total spin read from redata but not used
real(dp):: tt       ! Target temperature. Read in redata. Used in dynamics rout.
real(dp):: tstress(3,3) ! Total stress tensor, includin kinetic components
real(dp):: Uatm     ! Harris hartree electron energy, calculated in dhscf
real(dp):: Uscf     ! SCF hartree electron energy, calculated in dhscf
real(dp):: vcell(3,3) ! Unit cell "velocity" (time derivative)
real(dp):: virial   ! Virial = Sum(x*f)
real(dp):: vn       ! Velocity (time derivative) of the Nose thermostat
real(dp):: volcel   ! Function returning unit cell volume
real(dp):: volume   ! Unit cell volume
real(dp):: vpr      ! Velocity (time derivative) of the PR variables
real(dp):: wmix     ! Mixing weight for DM in SCF iteration
real(dp):: wmixkick ! Mixing weight for DM in special "kick" SCF steps
real(dp):: stressl(3,3) ! Local-mode part of the stress tensor
real(dp):: veclen   ! Length of a unit-cell vector

real(dp), allocatable :: &
  auxpul(:,:), &! Auxilialy matrices for Pulay mixing
  cfa(:,:),    &! Atomic forces orthogonalized to geometry constraints
  ebk(:,:,:),  &! Band energies
  efs(:),      &! Fermi energy for each spin (for fixed spin calculations)
  eo(:,:,:),   &! Hamiltonian eigenvalues
  fa(:,:),     &! Atomic forces
  fal(:,:),    &! Local-node part of atomic forces
  kpoint(:,:), &! Reciprocal vectors for BZ integration
  kweight(:),  &! Integration weights of BZ k-points
  polR(:,:),   &! Macroscopic polarization in unit-cell coordinates
  polxyz(:,:), &! Macroscopic polarization in cartesian coordinates
  qo(:,:,:),   &! Occupations of eigenstates
  qs(:),       &! Spin charges
  r2ij(:),     &! Interatomic distances squared
  va(:,:),     &! Atomic velocities
  xij(:,:),    &! Ineratomic vectors
  xijo(:,:),   &! Vectors between basis-orbital centers
  wfk(:,:)      ! Selected k points for eigenvalue printout

real(dp), pointer, save :: &
  bk(:,:),      &! Band k points
  Dold(:,:),    &! Density matrix of previous SCF step
  Dscf(:,:),    &! Density matrix of present SCF step
  Dscfsave(:,:),&! Auxiliary DM used in extrapol
  Eold(:,:),    &! Energy-density matrix of previous SCF step
  Escf(:,:),    &! Energy-density matrix of present SCF step
  H(:,:),       &! Hamiltonian matrix
  H0(:),        &! Non-SCF-dependent part of the Hamiltonian matrix
  kpol(:,:),    &! k points for the calculation of the polarization
  S(:),         &! Overlap matrix between basis orbitals
  wgthpol(:)     ! k-point integration weights for calculating the polarization

logical:: atmonly    ! Set up pseudoatom information only?
logical:: auxchanged ! Has auxiliary supercell changed?
logical:: bornz      ! Calculate Born polarization charges?
logical:: broyden_optim ! Use Broyden method to optimize geometry?
logical:: chebef     ! Compute the chemical potential in ordern?
logical:: default    ! Temporary used to pass default values in fdf reads
logical:: dminit     ! Initialize density matrix?
logical:: do_pdos    ! Compute the projected density of states?
logical:: dumpcharge ! Write electron density?
logical:: eggbox_block=.true. ! Read eggbox info from data file?
logical:: final=.false. ! Last geometry step?
logical:: first      ! Is this the first SCF step?
logical:: fixauxcell ! Keep the auxiliary supercell fixed?
logical:: fixspin    ! Keep the total spin fixed?
logical:: found      ! Found density matrix file?
logical:: foundxv    ! Found position-velocity file?
logical:: foundzm    ! Found Z-matrix file?
logical:: gamma      ! Gamma point (k=0) in BZ only?
logical:: genlogic   ! Generate local density of states?
logical:: harrisfun  ! Use Harris functional?
logical:: initdmaux  ! Re-initialize DM when auxiliary supercell changes?
logical:: inspn      ! Antiferro spin ordering in initdm?
logical:: itest      ! Temporary used to check if DM have just been read
logical:: last       ! Is this the last SCF step (done to calculate forces)?
logical:: lastst     ! Temporary to check if this is last geometry step
logical:: mix        ! Mix first SCF step? Used in broyden_mixing
logical:: mmix       ! Same as mix. Used in pulayx
logical:: muldeb     ! Write Mulliken polpulations at every SCF step?
logical:: naiveauxcell ! Use naive recipe for auxiliary supercell?
logical:: negl       ! Neglect hamiltonian matrix elements without overlap?
logical:: noeta      ! Use computed chemical potential instead of eta in ordern?
logical:: outlng     ! Long output?
logical:: overflow   ! Too short array sizes?
logical:: overflowed ! Have array sizes been increased?
logical:: pulfile    ! Use file to store Pulay info in pulayx?
logical:: relaxd=.false. ! Is geometry already relaxed?
logical:: require_energy_convergence ! to finish SCF iteration?
logical:: savdrh     ! Write file with diff. between SCF and atomic density?
logical:: savehs     ! Write file with Hamiltonian and overlap matrices?
logical:: savevh     ! Write file with Hartree electrostatic potential?
logical:: savevt     ! Write file with total effective potential?
logical:: savepsch   ! Write file with ionic (local pseudopotential) charge?
logical:: savetoch   ! Write file with total charge?
logical:: savrho     ! Write file with electron density?
logical:: struct_only ! Output initial structure only?
logical:: use_struct_file ! Read structural information from a special file?
logical:: usesavecg  ! Use continuation file for CG geometry relaxation?
logical:: usesavedm  ! Use cont. file for density matrix?
logical:: usesavedmloc ! Temporary to keep usesavedm value
logical:: usesavelwf ! Use continuation file for Wannier functions?
logical:: usesavexv  ! Use cont. file for atomic positions and velocities?
logical:: usesavezm  ! Use cont. file for Z-matrix?
logical:: varcel     ! Change unit cell during relaxation or dynamics?
logical:: writb      ! Write band eigenvalues?
logical:: writbk     ! Write k vectors of bands?
logical:: writedm    ! Write file with density of states?
logical:: writec     ! Write atomic coordinates at every geometry step?
logical:: writef     ! Write atomic forces at every geometry step?
logical:: writeig    ! Write eigenvalues?
logical:: writek     ! Write the k vectors of the BZ integration mesh?
logical:: writic     ! Write the initial atomic ccordinates?

integer, parameter :: llp5 = label_length+5
character(len=llp5):: fildrh   ! Filename for diff. between SCF and atomic dens.
character(len=llp5):: filepsch ! Filename for ionic (local pseudo) density
character(len=llp5):: filetoch ! Filename for total (ionic+electronic) dens.
character(len=llp5):: filevh   ! Filename for Hartree electrostatic potential
character(len=llp5):: filevt   ! Filename for total effective potential
character(len=llp5):: filrho   ! Filename for electron density
character(len=150) :: line     ! Data line read
character(len=10)  :: message  ! Not used
character(len=llp5):: paste    ! String concatenation function
character(len=150) :: sname    ! System name, used to initialise read
character(len=10)  :: shape    ! Geometry type (bulk,slab,wire,...)

type(parsed_line),pointer:: p ! Temporary to hold a digested line read

external:: automatic_cell ! Finds automatic unit cell
external:: bonds    ! Finds interatomic "bonds" (close atom pairs)
external:: cgvc     ! CG geometry relaxation including unit cell
external:: cgvc_zmatrix ! CG geometry relaxation using Z-matrix coordinates
external:: dhscf    ! Finds SCF-dependent part of the Hamiltonian matrix
external:: diagon   ! Solves the generalized eigenvalue problem
external:: dnaefs   ! Finds a correction to the neutral-atom energy&forces
external:: eggbox   ! Subtracts eggbox effect from energy and forces
external:: extrapol ! Extrapolates DM for new geometry
external:: fixed    ! Orthogonalizes atomic forces & stress to constraints
external:: initatom ! Initializes pseudoatom info
external:: iodm     ! I/O of density matrix
external:: ioeig    ! I/O of electron eigenenergies
external:: iofa     ! I/O of atomic forces
external:: iomd     ! I/O of atomic positions and velocities in MD
external:: iozm     ! I/O of Z-matrix
external:: kgrid    ! Finds k vectors for BZ integration
external:: kgridinit ! Initializes kgrid
external:: kinefsm  ! Finds kinetic-energy part of the Hamiltonian
external:: memory   ! Counts allocated memory sizes
external:: mulliken ! Finds Mulliken population charges
external:: naefs    ! Finds neutral-atom part of Hamiltonian, energy&forces
external:: neighb   ! Finds lists of neighbor atoms
external:: paste    ! Concatenates strings
external:: prversion ! Prints version of the program
external:: pulayx   ! Mixes density matrix using Pulay algorithm
external:: reinit   ! Initializes data read
external:: shaper   ! Finds geometry type
external:: spnvec   ! Finds spin vector
external:: timer    ! Measures and prints CPU time
external:: volcel   ! Finds unit-cell volume
external:: xijorb   ! Finds vectors between basis-orbital centers

#ifdef MPI
  integer :: MPIerror ! Return error code in MPI routines
  integer :: ntmp     ! Temporary to call an MPI routine
  real(dp):: buffer1  ! Temporary for the result of a globalize operation
  real(dp):: stresstmp(3,3) ! Temporary for globalized stress
  real(dp):: qtmp(4)  ! Temporary for globalized spin charge
  real(dp),allocatable:: fatmp(:,:) ! Temporary for globalized atomic forces
#endif

END MODULE m_siesta
