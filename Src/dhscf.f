C $Id: dhscf.f,v 1.29 1999/03/04 16:21:08 jose Exp $

      SUBROUTINE DHSCF( NSPIN, MAXORB, NORB, IAORB, IPHORB, INDXUO,
     .                  NA, ISA, XA, INDXUA, UCELL, MSCELL, G2MAX, NTM,
     .                  ILH, IFA, ISTR, IHMAT, FILRHO, FILDRH, FILEVH, 
     .                  FILEVT, MAXND, NUMD, LISTD, DSCF, DATM,
     .                  MAXNH, NUMH, LISTH, HMAT,
     .                  ENAATM, ENASCF, UATM, USCF, DUSCF, DUEXT, 
     .                  EXC, DXC,
     .                  DIPOL, FA, STRESS, IERROR )

C *******************************************************************
C FINDS SELFCONSISTENT-FIELD CONTRIBUTIONS TO HAMILTONIAN
C MATRIX ELEMENTS, TOTAL ENERGY AND ATOMIC FORCES.
C ALL DISTANCES AND ENERGIES IN ATOMIC UNITS (BOHRS AND RYDBERGS).
C CODED BY J.M.SOLER, AUGUST 1996. JULY 1997.
C    6  10        20        30        40        50        60        7072
C ************************* INPUT ***********************************
C INTEGER NSPIN         : NUMBER OF DIFFERENT SPIN POLARIZATIONS
C                         NSPIN=1 => UNPOLARIZED, NSPIN=2 => POLARIZED
C                         NSPIN=4 => Noncollinear spin
C INTEGER MAXORB        : SECOND DIMENSION OF DSCF AND HMAT
C INTEGER NORB          : Total number of basis orbitals in supercell
C INTEGER IAORB(NORB)   : ATOM TO WHICH EACH ORBITAL BELONGS
C INTEGER IPHORB(NORB)  : ORBITAL INDEX (WITHIN ATOM) OF EACH ORBITAL
C INTEGER INDEXUO       : Index of equivalent orbital in unit cell
C INTEGER NA            : Number of atoms in supercell
C INTEGER ISA(NA)       : Species idex of all atoms in supercell
C REAL*8  XA(3,NA)      : Atomic positions of all atoms in supercell
C INTEGER INDEXUA       : Index of equivalent atom in unit cell
C REAL*8  UCELL(3,3)    : UNIT CELL LATTICE VECTORS: UCELL(IXYZ,IVECT)
C INTEGER MSCELL(3,3)   : Supercell vectors in units of UCELL:
C                         SuperCell(IX,IV) = Sum_jv( UCELL(IX,JV) *
C                                                   MSCELL(JV,IV) )
C INTEGER ILH           : SWITCH WHICH FIXES WHETHER NUMH AND LISTH
C                         ARE INPUT OR OUTPUT ( ILH=0 => INPUT
C                                               ILH=1 => OUTPUT )
C INTEGER IFA           : SWITCH WHICH FIXES WHETHER THE SCF CONTRIB.
C                         TO ATOMIC FORCES IS CALCULATED AND ADDED TO FA
C INTEGER ISTR          : SWITCH WHICH FIXES WHETHER THE SCF CONTRIB.
C                         TO STRESS IS CALCULATED AND ADDED TO STRESS
C INTEGER IHMAT         : SWITCH WHICH FIXES WHETHER THE HMAT MATRIX
C                         ELEMENTS ARE CALCULATED OR NOT.
C CHARACTER*(*) FILRHO  : Name of file to saving the electron density
C                         (If blank => not saved)
C CHARACTER*(*) FILDRH  : Name of file to saving Delta Rho (Rho - Rho_atoms)
C                         (If blank => not saved)
C CHARACTER*(*) FILEVH  : Name of file to save electrostatic potential
C                         (If blank => not saved)
C CHARACTER*(*) FILEVT  : Name of file to save total potential
C                         (If blank => not saved)
C INTEGER MAXND             : FIRST DIMENSION OF LISTD AND DSCF
C INTEGER NUMD(NORB)        : NUMBER OF NONZERO DENSITY-MATRIX
C                             ELEMENTS FOR EACH MATRIX ROW
C INTEGER LISTD(MAXND,NORB) : NONZERO-DENSITY-MATRIX-ELEMENT COLUMN
C                             INDEXES FOR EACH MATRIX ROW
C REAL*8  DSCF(MAXND,MAXORB,NSPIN): SCF DENSITY-MATRIX ELEMENTS
C REAL*8  DATM(NORB)        : HARRIS DENSITY-MATRIX DIAGONAL ELEMENTS
C                             (ATOMIC OCCUPATION CHARGES OF ORBITALS)
C INTEGER MAXNH             : FIRST DIMENSION OF LISTH AND HMAT
C ****** INPUT OR OUTPUT (DEPENDING ON WHETHER MESH IS CALCULATED)**
C INTEGER NTM(3)            : NUMBER OF MESH DIVISIONS OF EACH CELL
C                             VECTOR, INCLUDING SUBGRID.
C ****** INPUT OR OUTPUT (DEPENDING ON ARGUMENT ILH) ***************
C INTEGER NUMH(NORB)        : NUMBER OF NONZERO HAMILTONIAN-MATRIX
C                             ELEMENTS FOR EACH MATRIX ROW
C INTEGER LISTH(MAXNH,NORB) : NONZERO-HAMILTONIAN-MATRIX-ELEMENT COLUMN
C                             INDEXES FOR EACH MATRIX ROW
C REAL*8  HMAT(MAXNH,MAXORB,NSPIN): Hamiltonian matrix in sparse form,
C                             to which are added the matrix elements
C                                 <ORB_I | DeltaV | ORB_J>, where
C                             DeltaV = Vna + Vxc(SCF) + 
C                                      Vhartree(RhoSCF-RhoHarris)
C                             This matrix is initialized only if ILH=1
C ************************* OUTPUT **********************************
C REAL*8  ENAATM : INTEGRAL OF VNA * RHOATM
C REAL*8  ENASCF : INTEGRAL OF VNA * RHOSCF
C REAL*8  UATM   : HARRIS HARTREE ELECTRON-INTERACTION ENERGY
C REAL*8  USCF   : SCF HARTREE ELECTRON-INTERACTION ENERGY
C REAL*8  DUSCF  : ELECTROSTATIC (HARTREE) ENERGY OF
C                    (RHOSCF - RHOATM) DENSITY
C REAL*8  DUEXT  : INTERACTION ENERGY WITH EXTERNAL ELECTRIC FIELD
C REAL*8  EXC    : SCF EXCHANCHE-CORRELATION ENERGY
C REAL*8  DXC    : SCF DOUBLE-COUNTING CORRECTION TO EXC
C                    DXC = INTEGRALOF( (EPSXC - VXC) * RHO )
C                    ALL ENERGIES IN RYDBERGS
C REAL*8  DIPOL(3): ELECTRIC DIPOLE (IN A.U.)
C                   ONLY WHEN THE SYSTEM IS A MOLECULE
C INTEGER IERROR : IERROR=0 => NO ERROR OCCURRED
C                  IERROR=1 => BAD INTERNAL DIMENSIONS. RECOMPILE
C *********************** INPUT AND OUTPUT **************************
C REAL*8  G2MAX      : EFFECTIVE PLANEWAVE CUTOFF IN RY. DETERMINES
C                        THE MESH DENSITY AND THE PRECISION OF INTEGRALS
C                          ON INPUT : VALUE REQUIRED
C                          ON OUTPUT: VALUE USED, WHICH MAY BE LARGER
C REAL*8  FA(3,NA)   : ATOMIC FORCES, TO WHICH THE SCF CONTRIBUTION
C                        IS ADDED BY THIS ROUTINE WHEN IFA=1.
C                        THE SCF CONTRIBUTION IS MINUS THE DERIVATIVE
C                        OF ( ENASCF - ENAATM + DUSCF + EXC ) WITH
C                        RESPECT TO ATOMIC POSITIONS, IN  RY/BOHR
C REAL*8 STRESS(3,3) : STRESS TENSOR, TO WHICH THE SCF CONTRIBUTION
C                      IS ADDED BY THIS ROUTINE WHEN IFA=1.
C                      THE SCF CONTRIBUTION IS MINUS THE DERIVATIVE OF
C                         ( ENASCF - ENAATM + DUSCF + EXC ) / VOLUME
C                      WITH RESPECT TO THE STRAIN TENSOR, IN RY.
C ******************* USER-DEFINED ROUTINES CALLED ******************
C The following routines must exist:
C --------------------------------
C REAL*8 FUNCTION RCUT(IS,IO)
C   Returns cutoff radius of orbitals and KB projectors.
C Input:
C     INTEGER IS : Species index
C     INTEGER IO : Orbital index
C --------------------------------
C SUBROUTINE PHIATM(IS,IO,R,PHI,GRPHI)
C    Finds values and gradients of:
C    a) basis orbitals (IO > 0)
C    b) KB proyectors  (IO < 0)
C    c) Local part of screened pseudopotential (IO = 0) ( b) and c) are
C       not required if MATEL is called only with IO > 0 )
C Input:
C   INTEGER IS   : Species index
C   INTEGER IO   : Orbital index
C   REAL*8  R(3) : Position with respect to atom
C Output:
C   REAL*8  PHI      : Value of orbital or KB projector at point R
C   REAL*8  GRPHI(3) : Gradient of PHI at point R
C --------------------------------
C REAL*8 FUNCTION RCORE(IS)
C   Returns cutoff radius of partial-core-correction charge density.
C   It must be zero if there is no partial-core-correction.
C Input:
C     INTEGER IS : Species index
C --------------------------------
C SUBROUTINE CHCORE(IS,R,CH,GRCH)
C   Finds partial-core-correction charge density and its gradient
C   This routine is called only if RCORE returns a nonzero cutoff 
C   radius for some species.
C Input:
C   INTEGER IS   : Species index
C   REAL*8  R(3) : Position with respect to atom
C Output:
C   REAL*8  CH      : Partial-core-correction charge density at point R
C   REAL*8  GRCH(3) : Gradient of CH at point R
C ************************ UNITS ************************************
C ENERGIES  IN RYDBERGS.
C DISTANCES IN BOHR.
C *******************************************************************

      IMPLICIT NONE

      INTEGER
     .  MAXND, MAXNH, MAXORB, NA, NORB, NSPIN,
     .  IAORB(NORB), IERROR, IFA, IHMAT, ILH, INDXUA(NA), INDXUO(NORB),
     .  IPHORB(NORB), ISA(NA), ISTR,
     .  LISTD(MAXND,NORB), LISTH(MAXNH,NORB),
     .  MSCELL(3,3), NTM(3), NUMD(NORB), NUMH(NORB)

      DOUBLE PRECISION
     .  DATM(MAXORB), DIPOL(3), DSCF(MAXND,MAXORB,NSPIN),
     .  DUEXT, DUSCF, DXC,
     .  ENAATM, ENASCF, EXC, FA(3,NA), G2MAX,
     .  HMAT(MAXNH,MAXORB,NSPIN),
     .  RCORE, RCUT, STRESS(3,3), UATM, UCELL(3,3), USCF, XA(3,NA)

      CHARACTER
     .  FILEVH*(*), FILEVT*(*), FILDRH*(*), FILRHO*(*)

      EXTERNAL
     .  CHCORE, PHIATM, RCORE, RCUT

C Fix number of sub-points (recomended NSM = 2) -------------------
      INTEGER NSM, NSP
      PARAMETER ( NSM = 2 )
      PARAMETER ( NSP = NSM * NSM * NSM )
C -----------------------------------------------------------------

C *******************************************************************
C Routines called internally:
C        CELLXC(...)    : Finds total exch-corr energy and potential
C        CHKDIM(...)    : Checks array dimensions
C        CHKGMX(...)    : Checks that mesh is fine enough for cutoff
C        CROSS(A,B,C)   : Finds the cross product of two vectors
C        DFSCF(...)     : Finds SCF contribution to atomic forces
C        DIGCEL(...)    : Finds a diagonal unit-cell/supercell 
C        DIPOLE(...)    : Finds electric dipole moment
C REAL*8 DISMIN(CELL,X) : Finds the minimum distance between point X
C                         and parallelepiped of CELL vectors
C REAL*8 DOT(V1,V2,N)   : Finds the dot product of vectors V1 and V2
C                         of length N.
C        EFIELD(...)    : Adds potential of an external electric field
C        IORHO(...)     : Saves electron density on a file
C        NFFT(N)        : Changes N into the next appropriate integer
C                         for the FFT routine
C        POISON(...)    : Solves Poisson equation
C        PRMEM(...)     : Prints array sizes
C        RECLAT(C,RC,I) : Finds reciprocal lattice vectors
C                         multiplied by 2*pi (I=1) or not (I=0).
C        REORD(...)     : Reorders electron density and potential arrays
C        RHOODA(...)    : Finds Harris electron density in the mesh
C        RHOOFD(...)    : Finds SCF electron density in the mesh
C        SHAPER(...)    : Finds the topology (shape) of the system
C        TIMER(...)     : Finds CPU times
C        TRANSP(...)    : Finds transpose of matrix PHI (basis
C                         orbitals ordered by mesh points
C        VLIST(...)     : Finds list of nonzero matrix elements
C                         of SCF potential
C        VMAT(...)      : Finds matrix elements of SCF potential
C REAL*8 VOLCEL( CELL ) : Returns volume of unit cell
C *******************************************************************
C Internal variables and arrays:
C REAL*8  AUX(*)        : General-purpose auxiliary-memory array
C LOGICAL BADDIM        : Bad dimensions?
C REAL*8  BCELL(3,3)    : Bulk lattice vectors
C REAL*8  CELL(3,3)     : Auxiliary lattice vectors (same as UCELL)
C REAL*8  CMESH(3,3)    : Mesh-cell vectors
C REAL*8  CONST         : Auxiliary variable (constant within a loop)
C REAL*8  DEC           : Auxiliary variable to call CELLXC
C REAL*8  DEX           : Auxiliary variable to call CELLXC
C REAL*8  DFA(3)        : Auxiliary variable to call DFSCF
C REAL*8  DFAATM(3)     : Auxiliary variable to call DFSCF
C REAL*8  DFASCF(3)     : Auxiliary variable to call DFSCF
C REAL*8  DSTRAT(3,3)   : Auxiliary variable to call DFSCF
C REAL*8  DSTRSC(3,3)   : Auxiliary variable to call DFSCF
C REAL*8  DVOL          : Mesh-cell volume
C REAL*8  DX(3)         : Vector from atom to mesh sub-point
C REAL*8  DXP(3)        : Vector from atom to mesh point
C REAL*8  DXSP(3)       : Vector from atom to mesh subpoint
C REAL*8  DXA(3,NA)     : Atom position within mesh-cell
C REAL*8  EC            : Correlation energy
C INTEGER ENDPH(0:NORB) : Last position occupied by orbitals in
C                         arrays PHI and LISTPH
C INTEGER ENDPHT(0:NMP) : Last position occupied by points in LSTPHT
C REAL*8  EX            : Exchange energy
C REAL*8  FIELD(3)      : External electric field
C LOGICAL FRSTME        : First time this routine is called?
C LOGICAL FOUND         : File found?
C REAL*8  G2MESH         : Effective planewave cutoff of mesh used
C REAL*8  GRPHI(3,NSP,MOP) : Gradient of basis orbitals at mesh points
C REAL*8  GRRHO(3)      : Gradient of density
C REAL*8  GRVA(3)       : Gradient of neutral-atom potential
C INTEGER I             : General-porpose index
C INTEGER I1,I2,I3      : Mesh indexes in each mesh direction
C INTEGER IA            : Atom index
C LOGICAL IBM           : Using IBM-ESSL library?
C INTEGER IDOP(MOP)     : Extended-mesh-index displacement of points
C                         within a sphere of radius RMAX
C INTEGER INDEXP(NEP)   : Translation from extended to normal mesh index
C INTEGER INDPHT(NTOP)  : Translation from point-ordered to
C                         orbital-ordered basis orbital storage.
C INTEGER IO            : Orbital index
C INTEGER IP            : Point index
C INTEGER IP0           : Point index of some origin
C INTEGER IPA(NA)       : Mesh cell in which atom is
C INTEGER IPHI          : Orbital index
C INTEGER IPHTY(NTY)    : Orbital index of each orbital type
C INTEGER IOP           : Index of Orbital Points
C INTEGER IS            : Species index
C LOGICAL ISDIAG        : Is supercell diagonal?
C LOGICAL ISEFLD        : Is there an external electric field?
C INTEGER ISP           : Sub-Point index
C INTEGER ISPIN         : Spin index
C INTEGER ISTY(NTY)     : Species index of each orbital type
C INTEGER ITY           : Orbital-type index
C INTEGER IUA           : Index of atom in the unit cell
C INTEGER IX(3)         : Mesh indexes in each mesh direction
C INTEGER J             : General-porpose index
C INTEGER J1,J2,J3      : Mesh indexes in each mesh direction
C INTEGER JOP           : Points-of-an-orbital index
C REAL*8  K0(3)         : Zero-vector argument for routine CHKGMX
C REAL*8  LASTC(3,3)    : Cell vectors of last call
C REAL*8  LASTRA        : Maximum orbital range in last call
C REAL*8  LASTXA(3,NA)  : Atomic positions of last call
C INTEGER LISTPH(NTOP)  : List of non-zero orbital points
C INTEGER LSTPHT(NTOP)  : List of non-zero orbitals at point
C INTEGER MOP           : Maximum number of non-zero Orbital Points
C INTEGER N             : Orbital-point index
C INTEGER NAUX          : Required size of auxiliary array
C INTEGER NBCELL        : Number of independent bulk lattice vectors
C INTEGER NCELLS        : Number of unit cells in supercell
C INTEGER NE(3)         : Number of mesh-Extension intervals 
C                         in each direction
C INTEGER NEM(3)        : Extended-mesh divisions in each direction
C INTEGER NEP           : Number of extended-mesh points
C INTEGER NM(3)         : Number of Mesh divisions of each cell vector
C INTEGER NMP           : Number of mesh points in unit cell
C INTEGER NMSC(3)       : Mesh divisions of each supercell vector
C INTEGER NOTY(NTY)     : Number of orbitals of each type
C INTEGER NPCC          : Partial core corrections? (0=no, 1=yes)
C INTEGER NSC(3)        : Number of unit-cells in each supercell direct.
C INTEGER NSD           : Number of diagonal spin values (1 or 2)
C INTEGER NSM           : Number of mesh sub-divisions in each direction
C INTEGER NSP           : Number of sub-points of each mesh point
C INTEGER NTP           : Number of mesh Total Points in unit cell
C                         (including subpoints)
C INTEGER NTPFFT        : Size of arrays required in FFT routine
C INTEGER NEP           : Number of Extended Points
C INTEGER NOTY(NTY)     : Number of orbitals of each type
C INTEGER NTOP          : Total number of nonzero orbital points
C INTEGER NTY           : Number of orbital TYpes
C INTEGER NUA           : Number of atoms in unit cell
C INTEGER NZERO(3)      : Auxiliary array with zeros, to call EFIELD
C REAL*4  PHI(NSP,NTOP) : Basis orbitals in mesh points (sparse format)
C REAL*8  PHIP          : Auxiliary variable to call PHIATM
C REAL*8  PLDIST        : Distance between mesh planes
C REAL*8  R             : Distance between atom and mesh point
C REAL*8  RA            : Cutoff radius of neutral-atom potential
C REAL*8  RCEL(3,3)     : Reciprocal cell vectors (with 2*pi factor)
C REAL*8  RCMESH(3,3)   : Reciprocal mesh-cell vectors
C                         (WITHOUT 2*pi factor)
C REAL*8  RCTY(NTY)     : Radius of each orbital type
C REAL*4  RHOATM(NTP)   : Harris electron density
C REAL*4  RHOPCC(NTP)   : Partial-core-correction density for xc
C REAL*4  DRHO(NTP)     : Selfconsistent electron density difference
C REAL*8  RHOP          : Density at one point
C REAL*8  RHOTOT        : Total density at one point
C REAL*8  RMAX          : Maximum orbital radius
C REAL*8  R2O           : Square of an orbital radius
C REAL*8  R2SP(NSP)     : Distance to subpoints
C LOGICAL SAMESH        : Same mesh of last call?
C LOGICAL SAMEXA        : Same atomic positions of last call?
C REAL*8  SCELL(3,3)    : Supercell vectors
C CHARACTER SHAPE*10    : Name of system shape
C REAL*8  TINY          : A small constant
C REAL*8  UHARRS        : Hartree energy of Harris electron density
C REAL*8  VA            : Neutral-atom potential
C REAL*4  VAUX(NTP)     : Auxiliary potential array
C REAL*8  VECMOD        : Vector modulus
C REAL*4  VNA(NTP)      : Sum of neutral-atom potentials
C REAL*8  VOLUME        : Unit cell volume
C REAL*4  VSCF(NTP)     : Hartree potential of selfconsistent density
C REAL*8  VXC           : Exchange-correlation potential
C LOGICAL WITHIN        : Is a mesh point within orbital range?
C CHARACTER XCAUTH*10   : Initials of xc version authors
C CHARACTER XCFUNC*10   : Name of xc functional
C REAL*8  XDOP(3,MOP)   : Vector to mesh points within RMAX
C REAL*8  XDSP(3,NSP)   : Vector to mesh sub-points
C REAL*8  XGRPHI(3,3,NSP,MAXOP) : Outer prod. of position*Gradient(phi)
C REAL*8  X0(3)         : Center of molecule
C *******************************************************************
C Dimension parameters of internal variables:
C INTEGER MAXA      : MAXimum number of Atoms
C INTEGER MAXAUX    : MAXimum size of AUXiliary array
C INTEGER MAXEP     : MAXimum number of Extended-mesh Points
C INTEGER MAXMP     : MAXimum number of Mesh Points in unit cell
C INTEGER MAXO      : MAXimum number of Orbitals
C INTEGER MAXOP     : MAXimum number of non-zero Orbital Points
C INTEGER MAXPCC    : MAXimum Partial Core Corrections (0 or 1)
C INTEGER MAXSPN    : MAXimum number of different spin polarizations
C INTEGER MAXTOP    : MAXimum number of non-zero Total Orbital Points
C INTEGER MAXTP     : MAXimum number of Total mesh Points
C INTEGER MAXTY     : MAXimum number of orbital TYpes
C *******************************************************************

C Dimension parameters for internal variables ---------------------
      INCLUDE 'dhscf.h'
*     INTEGER MAXA, MAXAUX, MAXEP, MAXMP, MAXO, MAXOP,
*     INTEGER MAXPCC, MAXSPN, MAXTOP, MAXTP, MAXTY
*     PARAMETER ( MAXA   = 1 )
*     PARAMETER ( MAXAUX = 1 )
*     PARAMETER ( MAXEP  = 1 )
*     PARAMETER ( MAXMP  = 1 )
*     PARAMETER ( MAXO   = 1 )
*     PARAMETER ( MAXOP  = 1 )
*     PARAMETER ( MAXPCC = 1 )
*     PARAMETER ( MAXSPN = 1 )
*     PARAMETER ( MAXTOP = 1 )
*     PARAMETER ( MAXTP  = 1 )
*     PARAMETER ( MAXTY  = 1 )
C -----------------------------------------------------------------

C Derived dimension parameters ------------------------------------
      INTEGER MINTY, MAXTYP
      PARAMETER ( MINTY  = 100 )
      PARAMETER ( MAXTYP = MINTY + MAXTY )
C -----------------------------------------------------------------

C Internal variable types and dimensions --------------------------
      INTEGER
     .  ENDPH(0:MAXO), ENDPHT(0:MAXMP),
     .  I, I1, I2, I3,
     .  IA, IDOP(MAXOP), INDEXP(MAXEP),
     .  INDPHT(MAXTOP), IO, IP, IP0, IPA(MAXA), IPHI, IPHTY(MAXTYP),
     .  IOP, IS, ISP, ISPIN, ISTY(MAXTYP), ITY, IUA, IX(3),
     .  J, J1, J2, J3, JOP,
     .  LISTPH(MAXTOP), LSTPHT(MAXTOP), MOP, 
     .  N, NAUX, NBCELL, NCELLS, NE(3), NEM(3), NEP, 
     .  NM(3), NMP, NMSC(3), NOTY(MAXTYP), NPCC, NSC(3), NSD, 
     .  NTOP, NTP, NTPFFT, NTY, NUA, NZERO(3)

      REAL
     .  PHI(NSP,MAXTOP), DRHO(MAXTP,MAXSPN), RHOATM(MAXTP),
     .  RHOPCC(MAXTP*MAXPCC+1), 
     .  VAUX(MAXTP), VNA(MAXTP), VSCF(MAXTP,MAXSPN)

      DOUBLE PRECISION
     .  AUX(MAXAUX), BCELL(3,3), B1XB2(3),
     .  CELL(3,3), CMESH(3,3), CONST,
     .  DEC, DEX, DFA(3), DFAATM(3), DFASCF(3),
     .  DISMIN, DOT, DSTRAT(3,3), DSTRES(3,3), DSTRSC(3,3),
     .  DVOL, DX(3), DXP(3), DXA(3,MAXA), DXSP(3,NSP),
     .  EC, EX, FIELD(3),
     .  G2MESH, GRPHI(3,NSP,MAXOP), GRRHO(3), GRVA(3),
     .  K0(3), LASTC(3,3), LASTRA, LASTXA(3,MAXA),
     .  PLDIST, PHIP,
*    .  QATM, QSCF,
     .  R, RA, RCELL(3,3), RCMESH(3,3), RCTY(MAXTYP),
     .  RHOP, RHOTOT, RMAX, R2O, R2SP(NSP), SCELL(3,3),
     .  TINY, UHARRS, VA, VECMOD, VOLUME, VOLCEL, VXC, 
     .  X0(3), XDOP(3,MAXOP), XDSP(3,NSP),
     .  XGRPHI(3,3,NSP,MAXOP)

      LOGICAL
     .  BADDIM, FOUND, FRSTME, ISDIAG, ISEFLD, SAMESH, SAMEXA, WITHIN

      CHARACTER
     .  SHAPE*10, XCAUTH*10, XCFUNC*10

      EXTERNAL
     .  CELLXC, CHKDIM, CHKGMX, CROSS,
     .  DFSCF, DIGCEL, DIPOLE, DISMIN, DOT,  
     .  EFIELD, IORHO, NFFT, POISON, PRMEM,
     .  RECLAT, REORD, RHOODA, RHOOFD,
     .  SHAPER, TIMER, TRANSP, VLIST, VMAT, VOLCEL

      SAVE
     .  CELL, CMESH, DVOL, DXA, ENDPH, ENDPHT, FIELD, FRSTME,
     .  G2MESH, K0,
     .  IDOP, INDEXP, INDPHT, IPA, ISDIAG, ISEFLD,
     .  LASTC, LASTRA, LASTXA, LISTPH, LSTPHT, MOP,
     .  NAUX, NE, NEM, NEP, NM, NMP, NMSC,
     .  NPCC, NSC, NTP, NTPFFT, NUA, NZERO,
     .  PHI, RCELL, RCMESH, RHOATM, RHOPCC, RMAX,
     .  SAMESH, SAMEXA, SCELL, SHAPE,
     .  TINY, UHARRS, VNA, VOLUME, XCAUTH, XCFUNC, XDOP, XDSP

      INCLUDE 'fdf/fdfdefs.h'

      DATA
     .  K0    / 3*0.D0 /,
     .  NZERO / 3*0 /,
     .  TINY  / 1.D-12 /,
     .  FRSTME /.TRUE./

C     Two unprobable numbers:
      DATA LASTC(1,1)  / 0.743978657912656D50 /
      DATA LASTXA(1,1) / 0.987654321273567D50 /
C -----------------------------------------------------------------

C Start time counter ----------------------------------------------
      CALL TIMER( 'DHSCF', 1 )
C -----------------------------------------------------------------

C Print array memory ----------------------------------------------
      IF (FRSTME) THEN
        CALL PRMEM( 0, 'DHSCF', 'AUX',    'D', MAXAUX        )
        CALL PRMEM( 0, 'DHSCF', 'DXA',    'D', 3*MAXA        )
        CALL PRMEM( 0, 'DHSCF', 'ENDPH',  'I', MAXO+1        )
        CALL PRMEM( 0, 'DHSCF', 'ENDPHT', 'I', MAXMP+1       )
        CALL PRMEM( 0, 'DHSCF', 'GRPHI',  'D', 3*NSP*MAXOP   )
        CALL PRMEM( 0, 'DHSCF', 'IDOP',   'I', MAXOP         )
        CALL PRMEM( 0, 'DHSCF', 'INDEXP', 'I', MAXEP         )
        CALL PRMEM( 0, 'DHSCF', 'INDPHT', 'I', MAXTOP        )
        CALL PRMEM( 0, 'DHSCF', 'IPA',    'I', MAXA          )
        CALL PRMEM( 0, 'DHSCF', 'IPHTY',  'I', MAXTYP        )
        CALL PRMEM( 0, 'DHSCF', 'ISTY',   'I', MAXTYP        )
        CALL PRMEM( 0, 'DHSCF', 'LASTXA', 'D', 3*MAXA        )
        CALL PRMEM( 0, 'DHSCF', 'LISTPH', 'I', MAXTOP        )
        CALL PRMEM( 0, 'DHSCF', 'LSTPHT', 'I', MAXTOP        )
        CALL PRMEM( 0, 'DHSCF', 'NOTY',   'I', MAXTYP        )
        CALL PRMEM( 0, 'DHSCF', 'PHI',    'R', NSP*MAXTOP    )
        CALL PRMEM( 0, 'DHSCF', 'RCTY',   'D', MAXTYP        )
        CALL PRMEM( 0, 'DHSCF', 'DRHO',   'R', MAXTP*NSPIN   )
        CALL PRMEM( 0, 'DHSCF', 'RHOATM', 'R', MAXTP         )
        CALL PRMEM( 0, 'DHSCF', 'RHOPCC', 'R', MAXTP*NPCC    )
        CALL PRMEM( 0, 'DHSCF', 'VAUX',   'R', MAXTP         )
        CALL PRMEM( 0, 'DHSCF', 'VNA',    'R', MAXTP         )
        CALL PRMEM( 0, 'DHSCF', 'VSCF',   'R', MAXTP*NSPIN   )
        CALL PRMEM( 0, 'DHSCF', 'XDOP',   'D', 3*MAXOP       )
        CALL PRMEM( 0, 'DHSCF', 'XGRPHI', 'D', 3*3*NSP*MAXOP )
        CALL PRMEM( 0, 'DHSCF', ' ',      ' ', 0             )
      ENDIF
C -----------------------------------------------------------------

C Initialize IERROR and 'bad dimensions' variable -----------------
      IERROR = 0
      BADDIM = .FALSE.
      IF ( NORB .GT. MAXO ) BADDIM = .TRUE.
C -----------------------------------------------------------------

C Find some data from the fdf input file -----------------
      IF (FRSTME) THEN
        XCFUNC = FDF_STRING('xc.functional','LDA')
        XCAUTH = FDF_STRING('xc.authors','PZ')
      ENDIF
C -----------------------------------------------------------------

C Find if mesh has to be changed ----------------------------------
      SAMESH = .TRUE.
      DO 20 I = 1,3
        DO 10 J = 1,3
          IF ( UCELL(J,I) .NE. LASTC(J,I) ) SAMESH = .FALSE.
          LASTC(J,I) = UCELL(J,I)
   10   CONTINUE
   20 CONTINUE
      IF ( G2MAX .GT. G2MESH * (1.D0 + TINY) ) SAMESH = .FALSE.
      RMAX = 0.D0
      DO 30 IO = 1,NORB
        IA = IAORB(IO)
        IPHI = IPHORB(IO)
        IS = ISA(IA)
        RMAX = MAX( RMAX, RCUT(IS,IPHI) )
   30 CONTINUE
      IF (RMAX .NE. LASTRA) SAMESH = .FALSE.
      LASTRA = RMAX
C ------------------------------------------------------------------

C Find if atoms have moved -----------------------------------------
      IF (NA .LE. MAXA) THEN
        SAMEXA = .TRUE.
        DO 50 IA = 1,NA
          DO 40 I = 1,3
            IF ( XA(I,IA) .NE. LASTXA(I,IA) ) SAMEXA = .FALSE.
            LASTXA(I,IA) = XA(I,IA)
   40     CONTINUE
   50   CONTINUE
      ELSE
        BADDIM = .TRUE.
        SAMEXA = .FALSE.
      ENDIF
C ------------------------------------------------------------------

C Mesh initialization ----------------------------------------------
      IF ( .NOT.SAMESH ) THEN

C       Start time counter for mesh initialization
*       CALL TIMER( 'DHSCF1', 1 )

C       Find diagonal unit cell and supercell
        CALL DIGCEL( UCELL, MSCELL, CELL, SCELL, NSC, ISDIAG )
        IF (.NOT.ISDIAG)
     .    WRITE(6,'(/,A,3(/,A,3F12.6,A,I6))')
     .      'DHSCF: WARNING: New shape of unit cell and supercell:',
     .     ('DHSCF:',(CELL(I,J),I=1,3),'   x',NSC(J),J=1,3)
        NCELLS = NSC(1) * NSC(2) * NSC(3)
        NUA = NA / NCELLS

C       Find shape of the system
        CALL SHAPER( CELL, NUA, ISA, XA, SHAPE, NBCELL, BCELL )

C       Find reciprocal cell vectors (multiplied by 2*pi)
        CALL RECLAT( CELL, RCELL, 1 )

C       Find number of mesh intervals for each cell vector.
C           Loop over cell vectors
        DO 60 I = 1,3
C             VECMOD is the legth of a reciprocal cell lattice vector
C             DOT makes a dot product of two vectors.
C             The reciprocal vectors of the mesh unit cell (CELL/NTM)
C             are RCELL*NTM, and must be larger than 2*GMAX
          VECMOD = SQRT( DOT( RCELL(1,I), RCELL(1,I), 3 ) )
          NTM(I) = 2 * SQRT(G2MAX) / VECMOD + 1
C             NFFT selects appropriate number of points for fft
   55     CALL NFFT( NTM(I) )
C             Require that NTM(I) to be a multiple of NSM
          IF ( MOD( NTM(I), NSM ) .NE. 0 ) THEN
            NTM(I) = NTM(I) + 1
            GOTO 55
          ENDIF
          NM(I) = NTM(I) / NSM
          NMSC(I) = NM(I) * NSC(I)
   60   CONTINUE
        WRITE(6,'(/,A,3(I6,A),I12)') 'DHSCF: MESH =',
     .    NTM(1),' x',NTM(2),' x',NTM(3),' =', NTM(1)*NTM(2)*NTM(3)

C       Find number of mesh points in unit cell.
        NMP = NM(1) * NM(2) * NM(3)
        NTP = NTM(1) * NTM(2) * NTM(3)
        NTPFFT = NTP + 2 * NTM(2) * NTM(3)
        IF (NMP    .GT. MAXMP) BADDIM = .TRUE.
        IF (NTPFFT .GT. MAXTP) BADDIM = .TRUE.

C       Find required size of auxiliary array REAL*8 AUX(NAUX)
C          Real*8  AUX(?) required in POISON
C          Integer AUX(NMP) required in TRANSP
C          Real*8  AUX(NORB) required in DFSCF
C          Real*4  AUX(NM1*NM2*NSP) required in REORD
        NAUX = -1
        CALL POISON( CELL, NTM(1), NTM(2), NTM(3), DRHO,
     .               DUSCF, VSCF, DSTRES, NAUX, AUX )
        NAUX = MAX( NAUX, (NMP+1) / 2 )
        NAUX = MAX( NAUX, NORB )
        NAUX = MAX( NAUX, (NM(1)*NM(2)*NSP+1) / 2 )
C       Next line is specific for fft.cray
*       NAUX = MAX( NAUX, 2 * ((NTM(1)+1)*NTM(2)*NTM(3)) )
        IF (NAUX .GT. MAXAUX) BADDIM = .TRUE.

C       Find volume of unit cell and of mesh cell
        VOLUME = VOLCEL( CELL )
        DVOL = VOLUME / NTP

C       Find effective cutoff
        G2MESH = 1.D6
        CALL CHKGMX( K0, RCELL, NTM, G2MESH )
        WRITE(6,'(A,2F10.3,A)')
     .   'DHSCF: Mesh cutoff (required, used) =', G2MAX, G2MESH, ' Ry'
        G2MAX = G2MESH

C       Find mesh-cell vectors
        DO 80 I = 1,3
          DO 70 J= 1,3
            CMESH(J,I) = CELL(J,I) / NM(I)
   70     CONTINUE
   80   CONTINUE

C       Find reciprocal mesh-cell vectors (not multiplied by 2*pi)
        CALL RECLAT( CMESH, RCMESH, 0 )

C       Find number of extended-mesh intervals for each cell vector.
C          Loop over mesh directions
        DO 90 I = 1,3
C            PLDIST is the distance between mesh planes
          PLDIST = 1.D0 / SQRT( DOT( RCMESH(1,I), RCMESH(1,I), 3 ) )
C            Find number of planes spanned by RMAX
          NE(I) = RMAX / PLDIST + DBLE( NSM-1 ) / DBLE( NSM )
          NE(I) = MAX( NE(I), NSM-2 )
C            Add NE(I) points to the left and NE(I)+1 points to the
C            right, to cover the spilling RMAX from an atom at any
C            possible place within the unit cell.
          NEM(I) = NMSC(I) + 2 * NE(I) + 1
   90   CONTINUE

C       Some printout for debugging
*       WRITE(6,'(A,3I6)') 'DHSCF: NM   =', NM 
*       WRITE(6,'(A,3I6)') 'DHSCF: NMSC =', NMSC 
*       WRITE(6,'(A,3I6)') 'DHSCF: NE   =', NE 
*       WRITE(6,'(A,3I6)') 'DHSCF: NEM  =', NEM 

C       Find total number of extended-mesh points.
        NEP = NEM(1) * NEM(2) * NEM(3)

C       Find relationship between extended and unit-cell mesh points
        IF (NEP .LE. MAXEP ) THEN
C         Loop over extended-mesh points
          DO 100 I3 = 0, NEM(3)-1
          DO 100 I2 = 0, NEM(2)-1
          DO 100 I1 = 0, NEM(1)-1
C           Find ext-mesh indexes in range [-NE(I),NM(I)+NE(I)]
            J1 = I1 - NE(1)
            J2 = I2 - NE(2)
            J3 = I3 - NE(3)
C           Find normal-mesh indexes in range [0,NM(I)]
C              1000*NM(I) is added to avoid negative numbers
C              in the argument of MOD
            J1 = MOD( J1 + 1000 * NMSC(1), NMSC(1) )
            J2 = MOD( J2 + 1000 * NMSC(2), NMSC(2) )
            J3 = MOD( J3 + 1000 * NMSC(3), NMSC(3) )
C              I = combined extended-mesh index.
            I = 1 + I1 + NEM(1) * I2 + NEM(1) * NEM(2) * I3
            IF (J1.LT.NM(1) .AND. J2.LT.NM(2) .AND. J3.LT.NM(3)) THEN
C               INDEXP(I) is the equivalent point within the unit cell
              INDEXP(I) = 1 + J1 + NM(1) * J2 + NM(1) * NM(2) * J3
            ELSE
              INDEXP(I) = -1
            ENDIF
  100     CONTINUE
        ELSE
          BADDIM = .TRUE.
        ENDIF

C       Find sub-points
        ISP = 0
        DO 140 I3 = 0, NSM-1
        DO 130 I2 = 0, NSM-1
        DO 120 I1 = 0, NSM-1
          ISP = ISP + 1
          DO 110 I = 1,3
            XDSP(I,ISP) = ( CMESH(I,1) * I1 +
     .                      CMESH(I,2) * I2 +
     .                      CMESH(I,3) * I3 ) / NSM
  110     CONTINUE
  120   CONTINUE
  130   CONTINUE
  140   CONTINUE

C       Find number of orbital types, number of orbitals of each type,
C       and cutoff radius of each type.
        NTY = 0
	DO ITY = 1,MAXTYP
	  NOTY(ITY) = 0
	ENDDO
        DO 160 IO = 1,NORB
          IA = IAORB(IO)
          IS = ISA(IA)
          IPHI = IPHORB(IO)
          DO 150 ITY = 1,NTY
            IF (IS.EQ.ISTY(ITY) .AND. IPHTY(ITY).EQ.IPHI) THEN
              NOTY(ITY) = NOTY(ITY) + 1
              GOTO 160
            ENDIF
  150     CONTINUE
          NTY = NTY + 1
          IF (NTY .GT. MAXTYP) STOP 'DHSCF: DIMENSION MAXTYP TOO SMALL'
          ISTY(NTY) = IS
          IPHTY(NTY) = IPHI
          NOTY(NTY) = 1
          RCTY(NTY) = RCUT( IS, IPHI )
  160   CONTINUE

C       Find points within RMAX (orbital points)
        MOP = 0
        NTOP = 0
C       Loop over possible points within RMAX
        DO 260 I3 = -NE(3), NE(3)+1
        DO 250 I2 = -NE(2), NE(2)+1
        DO 240 I1 = -NE(1), NE(1)+1
C         Find point coordinates
          DO 190 I = 1,3
            DXP(I) = CMESH(I,1) * I1 +
     .               CMESH(I,2) * I2 +
     .               CMESH(I,3) * I3 
  190     CONTINUE
C         Loop over sub-points
          WITHIN = .FALSE.
          DO 210 ISP = 1,NSP
C           Find point coordinates
            DO 200 I = 1,3
              DX(I) = DXP(I) + XDSP(I,ISP) 
  200       CONTINUE
C           Find distance from point to mesh cell
            R = DISMIN( CMESH, DX )
            IF ( R .LT. RMAX ) WITHIN = .TRUE.
  210     CONTINUE
          IF ( WITHIN ) THEN
C              MOP is the number of mesh points within RMAX
            MOP = MOP + 1
            DO 220 ITY = 1,NTY
              IF ( R .LT. RCTY(ITY) ) NTOP = NTOP + NOTY(ITY)
  220       CONTINUE
            IF ( MOP .LE. MAXOP ) THEN
C             Store index-distance and vector-distance to point.
              IDOP(MOP) = I1 + NEM(1) * I2 + NEM(1) * NEM(2) * I3
              DO 230 I = 1,3
                XDOP(I,MOP) = DXP(I)
  230         CONTINUE
            ELSE
              BADDIM = .TRUE.
            ENDIF
          ENDIF
  240   CONTINUE
  250   CONTINUE
  260   CONTINUE
        NTOP = NTOP / NCELLS
        IF (NTOP .GT. MAXTOP) BADDIM = .TRUE.

C       Find if there are partial-core-corrections for xc
        NPCC = 0
        DO 265 IA = 1,NA
          IF (RCORE(ISA(IA)) .GT. TINY) NPCC = 1
  265   CONTINUE
        IF (NPCC .LT. MAXPCC) BADDIM = .TRUE.

        IF (NSPIN .GT. MAXSPN) BADDIM = .TRUE.

C       Print good dimensions
        IF (BADDIM) THEN
          OPEN( 1, FILE='dhscf.h', STATUS='UNKNOWN' )
          WRITE(1,'(12(A,/))')
     .     'C DHSCF: Dimension parameters for DHSCF',
     .     'C MAXA   : MAXimum number of Atoms',
     .     'C MAXAUX : MAXimum size of AUXiliary array',
     .     'C MAXEP  : MAXimum number of Extended-mesh Points',
     .     'C MAXMP  : MAXimum number of Mesh Points in unit cell',
     .     'C MAXO   : MAXimum number of Orbitals',
     .     'C MAXOP  : MAXimum number of non-zero Orbital Points',
     .     'C MAXPCC : MAXimum Partial Core Corrections (0 or 1)',
     .     'C MAXSPN : MAXimum number of different spin polarizations',
     .     'C MAXTOP : MAXimum number of non-zero Total Orbital Points',
     .     'C MAXTP  : MAXimum number of Total mesh Points',
     .     'C MAXTY  : MAXimum number of orbital TYpes'
          WRITE(1,'(6X,A)')
     .      'INTEGER MAXA, MAXAUX, MAXEP, MAXMP, MAXO, MAXOP',
     .      'INTEGER MAXPCC, MAXSPN, MAXTOP, MAXTP, MAXTY'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXAUX =', NAUX,   ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXA   =', NA,     ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXEP  =', NEP,    ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXMP  =', NMP,    ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXO   =', NORB,   ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXOP  =', MOP,    ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXPCC =', NPCC,   ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXSPN =', NSPIN,  ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXTOP =', NTOP,   ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXTP  =', NTPFFT, ' )'
          WRITE(1,'(6X,A,I10,A)') 'PARAMETER ( MAXTY  =', NTY,    ' )'
          WRITE(1,'(/,A,F10.3,A)')'C DHSCF: Mesh cutoff =', G2MAX,' Ry'
          CALL PRMEM( 1, 'DHSCF', 'AUX',    'D', NAUX         )
          CALL PRMEM( 1, 'DHSCF', 'DXA',    'D', 3*NA         )
          CALL PRMEM( 1, 'DHSCF', 'ENDPH',  'I', NORB+1       )
          CALL PRMEM( 1, 'DHSCF', 'ENDPHT', 'I', NMP+1        )
          CALL PRMEM( 1, 'DHSCF', 'GRPHI',  'D', 3*NSP*MOP    )
          CALL PRMEM( 1, 'DHSCF', 'IDOP',   'I', MOP          )
          CALL PRMEM( 1, 'DHSCF', 'INDEXP', 'I', NEP          )
          CALL PRMEM( 1, 'DHSCF', 'INDPHT', 'I', NTOP         )
          CALL PRMEM( 1, 'DHSCF', 'IPA',    'I', NA           )
          CALL PRMEM( 1, 'DHSCF', 'IPHTY',  'I', NTY          )
          CALL PRMEM( 1, 'DHSCF', 'ISTY',   'I', NTY          )
          CALL PRMEM( 1, 'DHSCF', 'LASTXA', 'D', 3*NA         )
          CALL PRMEM( 1, 'DHSCF', 'LISTPH', 'I', NTOP         )
          CALL PRMEM( 1, 'DHSCF', 'LSTPHT', 'I', NTOP         )
          CALL PRMEM( 1, 'DHSCF', 'NOTY',   'I', NTY          )
          CALL PRMEM( 1, 'DHSCF', 'PHI',    'R', NSP*NTOP     )
          CALL PRMEM( 1, 'DHSCF', 'RCTY',   'D', NTY          )
          CALL PRMEM( 1, 'DHSCF', 'DRHO',   'R', NTPFFT*NSPIN )
          CALL PRMEM( 1, 'DHSCF', 'RHOATM', 'R', NTPFFT       )
          CALL PRMEM( 1, 'DHSCF', 'RHOPCC', 'R', NTPFFT*NPCC  )
          CALL PRMEM( 1, 'DHSCF', 'VAUX',   'R', NTPFFT       )
          CALL PRMEM( 1, 'DHSCF', 'VNA',    'R', NTPFFT       )
          CALL PRMEM( 1, 'DHSCF', 'VSCF',   'R', NTPFFT*NSPIN )
          CALL PRMEM( 1, 'DHSCF', 'XDOP',   'D', 3*MOP        )
          CALL PRMEM( 1, 'DHSCF', 'XGRPHI', 'D', 3*3*NSP*MOP  )
          CALL PRMEM( 1, 'DHSCF', ' ',      ' ', 0            )
          CLOSE(1)
*         STOP 'DHSCF: BAD DIMENSIONS. RECOMPILE.'
          IERROR = 1
          GOTO 999
        ENDIF

C       Stop time counter for mesh initialization
*       CALL TIMER( 'DHSCF1', 2 )

      ENDIF
C End of mesh initialization ------------------------------------------

C Initialize atomic orbitals, density and potential -------------------
      IF (.NOT.SAMEXA) THEN

C       Start time counter for atomic initializations
        CALL TIMER( 'DHSCF2', 1 )

C       Find atomic positions relative to mesh
        DO 290 IA = 1,NA
C         Find index of extended-mesh cell in which atom is
          DO 270 I = 1,3
            DX(I) = DOT( XA(1,IA), RCMESH(1,I), 3)
            IX(I) = INT( DX(I) + 100000 ) - 100000
            DX(I) = DX(I) - IX(I)
            IX(I) = MOD( IX(I) + 1000*NMSC(I), NMSC(I) )
            IX(I) = IX(I) + NE(I)
  270     CONTINUE
          IPA(IA) = 1 + IX(1) + NEM(1) * IX(2) +
     .              NEM(1) * NEM(2) * IX(3)
C         Find atom position within mesh cell
          DO 280 I = 1,3
            DXA(I,IA) = CMESH(I,1) * DX(1) +
     .                  CMESH(I,2) * DX(2) +
     .                  CMESH(I,3) * DX(3)
  280     CONTINUE
  290   CONTINUE

C       Find partial-core-correction energy density
        IF (NPCC .EQ. 1) THEN
          DO 292 IP = 1,NTP
            RHOPCC(IP) = 0.D0
  292     CONTINUE

          DO 300 IA = 1,NA
            IS = ISA(IA)
            RA = RCORE( IS )
            IF (RA .GT. TINY) THEN

C             Loop over mesh points inside RMAX
              DO 298 IOP = 1,MOP
                IP0 = INDEXP( IPA(IA) + IDOP(IOP) )
                IF (IP0 .GT. 0) THEN

C                 Loop over sub-points
                  DO 296 ISP = 1,NSP
                    DO 294 I = 1,3
                      DX(I) = XDOP(I,IOP) + XDSP(I,ISP) - DXA(I,IA)
  294               CONTINUE
                    R = SQRT( DOT( DX, DX, 3 ) )
                    IF (R .LT. RA) THEN
                      IP = ISP + NSP * (IP0 - 1)
                      CALL CHCORE( IS, DX, RHOP, GRRHO )
                      RHOPCC(IP) = RHOPCC(IP) + RHOP
                    ENDIF
  296             CONTINUE

                ENDIF
  298         CONTINUE

            ENDIF
  300     CONTINUE
        ENDIF

C       Find neutral-atom potential
        DO 305 IP = 1,NTP
          VNA(IP) = 0.D0
  305   CONTINUE

        DO 340 IA = 1,NA
          IS = ISA(IA)
          RA = RCUT( IS, 0 )

C         Loop over mesh points inside RMAX
          DO 330 IOP = 1,MOP
            IP0 = INDEXP( IPA(IA) + IDOP(IOP) )
            IF (IP0 .GT. 0) THEN

C             Loop over sub-points
              DO 320 ISP = 1,NSP
                DO 310 I = 1,3
                  DX(I) = XDOP(I,IOP) + XDSP(I,ISP) - DXA(I,IA)
  310           CONTINUE
                R = SQRT( DOT( DX, DX, 3 ) )
                IF (R .LT. RA) THEN
                  IP = ISP + NSP * (IP0 - 1)
                  CALL PHIATM( IS, 0, DX, VA, GRVA )
                  VNA(IP) = VNA(IP) + VA
                ENDIF
  320         CONTINUE

            ENDIF
  330     CONTINUE
  340   CONTINUE

C       Find atomic orbitals at mesh points
C       Loop over orbitals
        N = 0
        ENDPH(0) = 0
        DO 400 IO = 1,NORB
          IA = IAORB(IO)
          IPHI = IPHORB(IO)
          IS = ISA(IA)
          R2O = RCUT(IS,IPHI)**2

C         Loop over mesh points inside RMAX
          DO 390 IOP = 1,MOP
            IP0 = INDEXP( IPA(IA) + IDOP(IOP) )
            IF (IP0 .GT. 0) THEN

C             Loop over sub-points to find if point is within range
              WITHIN = .FALSE.
              DO 360 ISP = 1,NSP
                DO 350 I = 1,3
                  DXSP(I,ISP) = XDOP(I,IOP) + XDSP(I,ISP) - DXA(I,IA)
  350           CONTINUE
                R2SP(ISP) = DXSP(1,ISP)**2 + DXSP(2,ISP)**2 +
     .                      DXSP(3,ISP)**2
                IF (R2SP(ISP) .LT. R2O) WITHIN = .TRUE.
  360         CONTINUE

C             If within range, add point to list of orbital points
              IF (WITHIN) THEN
                N = N + 1
                CALL CHKDIM( 'DHSCF', 'MAXTOP', MAXTOP, N, 1 )
                LISTPH(N) = IP0
C               Loop again over sub-points to calculate PHI
                DO 380 ISP = 1,NSP
                  IF (R2SP(ISP) .LT. R2O) THEN
                    CALL PHIATM( IS, IPHI, DXSP(1,ISP),
     .                           PHIP, GRPHI )
                    PHI(ISP,N) = PHIP
                  ELSE
                    PHI(ISP,N) = 0.D0
                  ENDIF
  380           CONTINUE
              ENDIF

            ENDIF
  390     CONTINUE
          ENDPH(IO) = N
  400   CONTINUE

C       Find transpose of matrix PHI
        CALL TRANSP( NORB, ENDPH, LISTPH,
     .               NMP, ENDPHT, LSTPHT, INDPHT, 2*MAXAUX, AUX )

C       Find (ILH=1) or check (only if NSP=1) NUMH and LISTH
        IF (ILH.EQ.1 .OR. NSP.EQ.1)
     .    CALL VLIST( ILH, NORB,  ENDPH,  LISTPH,
     .                NMP,   ENDPHT, LSTPHT,
     .                MAXNH, NUMH,   LISTH, 2*MAXAUX, AUX )

C       Find Harris (sum of atomic) electron density
        CALL RHOODA( NORB, INDXUO, ENDPH, LISTPH, PHI,
     .               NSP, NSP, NMP, DATM, RHOATM )

C       Find Hartree energy of RHOATM, using VSCF as an auxiliary array
C       Reorder RHOATM into a sequential array in the total mesh
        CALL REORD( RHOATM, RHOATM, NM, NSM, +1, 2*MAXAUX, AUX )
C       Solve Poisson's equation
        NAUX = MAXAUX
        CALL POISON( CELL, NTM(1), NTM(2), NTM(3), RHOATM,
     .               UHARRS, VSCF, DSTRES, NAUX, AUX )
C       Reorder back RHOATM into mesh points and sub-points
        CALL REORD( RHOATM, RHOATM, NM, NSM, -1, 2*MAXAUX, AUX )

C       Stop time counter for atomic initializations
        CALL TIMER( 'DHSCF2', 2 )

      ENDIF
C End of atomic initializations ---------------------------------------

C Start time counter for SCF iteration part ---------------------------
      CALL TIMER( 'DHSCF3', 1 )
C ---------------------------------------------------------------------

C Initialize HMAT -----------------------------------------------------
      IF ( ILH .EQ. 1 ) THEN
        DO 415 ISPIN = 1,NSPIN
          DO 410 IO = 1,NORB
            DO 405 I = 1,NUMH(IO)
              HMAT(I,IO,ISPIN) = 0.D0
  405       CONTINUE
  410     CONTINUE
  415   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Find number of diagonal spin values ---------------------------------
      NSD = MIN( NSPIN, 2 )
C ---------------------------------------------------------------------

C Find SCF electron density at mesh points. Store it in array DRHO ----
      DO 420 ISPIN = 1,NSPIN
        CALL RHOOFD( NORB, INDXUO, ENDPH,  LISTPH, PHI, NSP, NSP,
     .               NMP,  ENDPHT, LSTPHT, INDPHT,
     .               MAXND, NUMD, LISTD, DSCF(1,1,ISPIN),
     .               DRHO(1,ISPIN), 2*MAXAUX, AUX )
  420 CONTINUE
C ---------------------------------------------------------------------

C Save electron density -----------------------------------------------
      IF (FILRHO .NE. ' ') THEN
        DO 422 ISPIN = 1,NSPIN
          CALL REORD( DRHO(1,ISPIN), DRHO(1,ISPIN),
     .                NM, NSM, +1, 2*MAXAUX, AUX )
  422   CONTINUE
        CALL IORHO( 'WRITE', FILRHO, CELL, NTM, MAXTP, NSPIN,
     .              DRHO, FOUND )
        DO 424 ISPIN = 1,NSPIN
          CALL REORD( DRHO(1,ISPIN), DRHO(1,ISPIN),
     .                NM, NSM, -1, 2*MAXAUX, AUX )
  424   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Find difference between selfconsistent and atomic densities ---------
      DO 430 ISPIN = 1,NSD
        DO 428 IP = 1,NTP
          DRHO(IP,ISPIN) = DRHO(IP,ISPIN) - RHOATM(IP) / NSD
  428   CONTINUE
  430 CONTINUE
C ---------------------------------------------------------------------

C Save electron density difference ------------------------------------
      IF (FILDRH .NE. ' ') THEN
        DO 432 ISPIN = 1,NSPIN
          CALL REORD( DRHO(1,ISPIN), DRHO(1,ISPIN),
     .                NM, NSM, +1, 2*MAXAUX, AUX )
  432   CONTINUE
        CALL IORHO( 'WRITE', FILDRH, CELL, NTM, MAXTP, NSPIN,
     .              DRHO, FOUND )
        DO 434 ISPIN = 1,NSPIN
          CALL REORD( DRHO(1,ISPIN), DRHO(1,ISPIN),
     .                NM, NSM, -1, 2*MAXAUX, AUX )
  434   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Find atomic and SCF total charges for debugging ---------------------
*     QATM = 0.D0
*     QSCF = 0.D0
*     DO 450 ISPIN = 1,NSD
*       QATM = 0.D0
*       QSCF = 0.D0
*       DO 440 IP = 1,NTP
*         QATM = QATM + DVOL * RHOATM(IP) / NSD
*         QSCF = QSCF + DVOL * ( RHOATM(IP) / NSD + DRHO(IP,ISPIN) )
* 440   CONTINUE
*       IF (NSD .EQ. 1) THEN
*         WRITE(6,'(A,2F12.6)') 'DHSCF: QATM,QSCF =', QATM, QSCF
*       ELSE
*         WRITE(6,'(A,I3,2F12.6)')
*    .       'DHSCF: ISPIN, QATM,QSCF =', ISPIN, QATM, QSCF
*       ENDIF
* 450 CONTINUE
C ---------------------------------------------------------------------

C Transform spin density into sum and difference ---------------------
      IF (NSD .EQ. 2) THEN
        I1 = 1
        I2 = 2
        DO 460 IP = 1,NTP
          RHOTOT = DRHO(IP,I1) + DRHO(IP,I2)
          DRHO(IP,I2) = DRHO(IP,I2) - DRHO(IP,I1)
          DRHO(IP,I1) = RHOTOT
  460   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Find electric dipole ------------------------------------------------
      IF (SHAPE .NE. 'bulk') THEN

C       Find center of system
        DO 463 I = 1,3
          X0(I) = 0.D0
          DO 462 IA = 1,NUA
            X0(I) = X0(I) + XA(I,IA) / NUA
  462     CONTINUE
  463   CONTINUE

C       Find dipole
        CALL REORD( DRHO, DRHO, NM, NSM, +1, 2*MAXAUX, AUX )
        CALL DIPOLE( CELL, NTM(1), NTM(2), NTM(3), DRHO, X0, DIPOL )
        CALL REORD( DRHO, DRHO, NM, NSM, -1, 2*MAXAUX, AUX )

C       Orthogonalize dipole to bulk directions
        IF (SHAPE .EQ. 'chain') THEN
          CONST = DOT( DIPOL, BCELL, 3 ) / DOT( BCELL, BCELL, 3 )
          DO 464 I = 1,3
            DIPOL(I) = DIPOL(I) - CONST * BCELL(I,1)
  464     CONTINUE
        ELSEIF (SHAPE .EQ. 'slab') THEN
          CALL CROSS( BCELL(1,1), BCELL(1,2), B1XB2 )
          CONST = DOT( DIPOL, B1XB2, 3 ) / DOT( B1XB2, B1XB2, 3 )
          DO 465 I = 1,3
            DIPOL(I) = CONST * B1XB2(I)
  465     CONTINUE
        ENDIF
      ENDIF
C ---------------------------------------------------------------------

C Find Hartree potential of DRHO = RHOSCF-RHOATM. Store it in array VSCF
C     Reorder DRHO into a sequential array in the total mesh
      CALL REORD( DRHO, DRHO, NM, NSM, +1, 2*MAXAUX, AUX )
C     Solve Poisson's equation
      DO 466 IP = 1,NTP
        VSCF(IP,1) = 0.D0
  466 ENDDO
      NAUX = MAXAUX
      CALL POISON( CELL, NTM(1), NTM(2), NTM(3), DRHO,
     .             DUSCF, VSCF, DSTRES, NAUX, AUX )
C     Reorder back DRHO and VSCF into mesh points and sub-points
      CALL REORD( DRHO, DRHO, NM, NSM, -1, 2*MAXAUX, AUX )
      CALL REORD( VSCF, VSCF, NM, NSM, -1, 2*MAXAUX, AUX )
C ---------------------------------------------------------------------

C Add contribution to stress from electrostatic energy of RHOSCF-RHOATM
      IF (ISTR .EQ. 1) THEN
        DO 470 J = 1,3
          DO 468 I = 1,3
            STRESS(I,J) = STRESS(I,J) + DSTRES(I,J)
  468     CONTINUE
  470   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Find electrostatic (Hartree) energy of full SCF electron density ---
      UATM = UHARRS
      USCF = UATM + DUSCF
      DO 475 IP = 1,NTP
        USCF = USCF + DVOL * VSCF(IP,1) * RHOATM(IP)
  475 CONTINUE
C ---------------------------------------------------------------------

C Add neutral-atom potential to VSCF ----------------------------------
      ENAATM = 0.D0
      ENASCF = 0.D0
      DO 480 IP = 1,NTP
        ENAATM = ENAATM + DVOL * VNA(IP) * RHOATM(IP)
        ENASCF = ENASCF + DVOL * VNA(IP) * ( RHOATM(IP) + DRHO(IP,1) )
        VSCF(IP,1) = VSCF(IP,1) + VNA(IP)
  480 CONTINUE
C ---------------------------------------------------------------------

C Add the potential of a (possible) external electric field -----------
      IF (FRSTME) THEN
        CALL EFIELD( CELL, NUA, ISA, XA, NZERO, VSCF, FIELD )
        ISEFLD = .FALSE.
        IF (SQRT(DOT(FIELD,FIELD,3)) .GT. TINY) ISEFLD = .TRUE.
      ENDIF
      IF (ISEFLD) THEN
        CALL REORD( VSCF, VSCF, NM, NSM, +1, 2*MAXAUX, AUX )
        CALL EFIELD( CELL, NUA, ISA, XA, NTM, VSCF, FIELD )
        CALL REORD( VSCF, VSCF, NM, NSM, -1, 2*MAXAUX, AUX )
        DUEXT = - DOT( FIELD, DIPOL, 3 )
      ENDIF
C ---------------------------------------------------------------------

C Save electrostatic potential ----------------------------------------
      IF (FILEVH .NE. ' ') THEN
        CALL REORD( VSCF, VSCF, NM, NSM, +1, 2*MAXAUX, AUX )
        CALL IORHO( 'WRITE', FILEVH, CELL, NTM, MAXTP, 1, VSCF, FOUND )
        CALL REORD( VSCF, VSCF, NM, NSM, -1, 2*MAXAUX, AUX )
      ENDIF
C ---------------------------------------------------------------------

C Add contribution to stress from the derivative of the Jacobian of ---
C r->r' (strained r) in the integral of VNA*(RHOSCF-RHOATM)
        IF (ISTR .EQ. 1) THEN
          DO 485 I = 1,3
            STRESS(I,I) = STRESS(I,I) + ( ENASCF - ENAATM ) / VOLUME
  485     CONTINUE
        ENDIF
C ---------------------------------------------------------------------

C Get back spin density from sum and difference ---------------------
      IF (NSD .EQ. 2) THEN
        I1 = 1
        I2 = 2
        DO 490 IP = 1,NTP
          RHOTOT = DRHO(IP,I1)
          DRHO(IP,I1) = 0.5D0 * ( RHOTOT - DRHO(IP,I2) )
          DRHO(IP,I2) = 0.5D0 * ( RHOTOT + DRHO(IP,I2) )
  490   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Find exchange-correlation energy and potential ----------------------
      DO 520 IP = 1,NTP
        VAUX(IP) = VSCF(IP,1)
  520 CONTINUE
      DO 535 ISPIN = 1,NSD
        DO 530 IP = 1,NTP
          DRHO(IP,ISPIN) = DRHO(IP,ISPIN) + RHOATM(IP) / NSD
          IF (NPCC .EQ. 1) 
     .      DRHO(IP,ISPIN) = DRHO(IP,ISPIN) + RHOPCC(IP) / NSD
  530   CONTINUE
  535 CONTINUE
      DO 540 ISPIN = 1,NSPIN
        CALL REORD(DRHO(1,ISPIN),DRHO(1,ISPIN),NM,NSM,+1,2*MAXAUX,AUX)
  540 CONTINUE
      CALL CELLXC( XCFUNC, XCAUTH, 0,
     .             CELL, NTM, NTM, MAXTP, 0, AUX, NSPIN, DRHO,
     .             EX, EC, DEX, DEC, VSCF, DSTRES, MAXAUX, AUX )
      EXC = EX + EC
      DXC = DEX + DEC
      DO 545 ISPIN = 1,NSPIN
        CALL REORD(DRHO(1,ISPIN),DRHO(1,ISPIN),NM,NSM,-1,2*MAXAUX,AUX)
        CALL REORD(VSCF(1,ISPIN),VSCF(1,ISPIN),NM,NSM,-1,2*MAXAUX,AUX)
  545 CONTINUE
      DO 560 ISPIN = 1,NSD
        DO 550 IP = 1,NTP
          DRHO(IP,ISPIN) = DRHO(IP,ISPIN) - RHOATM(IP) / NSD
          IF (NPCC .EQ. 1) 
     .      DRHO(IP,ISPIN) = DRHO(IP,ISPIN) - RHOPCC(IP) / NSD
          VSCF(IP,ISPIN) = VSCF(IP,ISPIN) + VAUX(IP)
  550   CONTINUE
  560 CONTINUE
C ---------------------------------------------------------------------

C Save total potential ----------------------------------------
      IF (FILEVT .NE. ' ') THEN
        DO 562 ISPIN = 1,NSPIN
          CALL REORD( VSCF(1,ISPIN), VSCF(1,ISPIN),
     .                NM, NSM, +1, 2*MAXAUX, AUX )
  562   CONTINUE
        CALL IORHO( 'WRITE', FILEVT, CELL, NTM, MAXTP, NSPIN,
     .              VSCF, FOUND )
        DO 564 ISPIN = 1,NSPIN
          CALL REORD( VSCF(1,ISPIN), VSCF(1,ISPIN),
     .                NM, NSM, -1, 2*MAXAUX, AUX )
  564   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Add contribution to stress from exchange-correlation energy ---------
      IF (ISTR .EQ. 1) THEN
        DO 574 J = 1,3
          DO 572 I = 1,3
            STRESS(I,J) = STRESS(I,J) + DSTRES(I,J)
  572     CONTINUE
  574   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Find SCF contribution to hamiltonian matrix elements ----------------
      IF (IHMAT .EQ. 1) THEN
        DO 580 ISPIN = 1,NSPIN
          CALL VMAT( NORB, INDXUO, ENDPH,  LISTPH, PHI, NSP, NSP,
     .               NMP,  ENDPHT, LSTPHT, INDPHT,
     .               VOLUME, VSCF(1,ISPIN),
     .               MAXNH, NUMH, LISTH, HMAT(1,1,ISPIN),
     .               2*MAXAUX, AUX )
  580   CONTINUE
      ENDIF
C ---------------------------------------------------------------------

C Stop time counter for SCF iteration part ----------------------------
      CALL TIMER( 'DHSCF3', 2 )
C ---------------------------------------------------------------------

C Find SCF contribution to atomic forces and/or stress ----------------
      IF (IFA.EQ.1 .OR. ISTR.EQ.1) THEN

C       Start time counter for force calculation part
        CALL TIMER( 'DHSCF4', 1 )

C       Transform spin density into sum and difference
        IF (NSD .EQ. 2) THEN
          I1 = 1
          I2 = 2
          DO 590 IP = 1,NTP
            RHOTOT = DRHO(IP,I1) + DRHO(IP,I2)
            DRHO(IP,I2) = DRHO(IP,I2) - DRHO(IP,I1)
            DRHO(IP,I1) = RHOTOT
  590     CONTINUE
        ENDIF

C       Find contribution of neutral-atom potential
        DO 650 IA = 1,NA
          IUA = INDXUA(IA)
          IS = ISA(IA)
          RA = RCUT( IS, 0 )
C         Loop over mesh points and subpoints inside RMAX
          DO 640 IOP = 1,MOP
            IP0 = INDEXP( IPA(IA) + IDOP(IOP) )
            IF (IP0 .GT. 0) THEN
              DO 630 ISP = 1,NSP
                DO 600 I = 1,3
                  DX(I) = XDOP(I,IOP) + XDSP(I,ISP) - DXA(I,IA)
  600           CONTINUE
                R = SQRT( DOT( DX, DX, 3 ) )
                IF (R .LT. RA) THEN
                  IP = ISP + NSP * (IP0 - 1)
                  CALL PHIATM( IS, 0, DX, VA, GRVA )
                  DO 620 I = 1,3
                    DFA(I) = DVOL * GRVA(I) * DRHO(IP,1)
                    IF (IFA .EQ. 1) FA(I,IUA) = FA(I,IUA) + DFA(I)
                    IF (ISTR .EQ. 1) THEN
*                     STRESS(I,I) = STRESS(I,I) +
*    .                              DVOL * VA * DRHO(IP,1) / VOLUME
                      DO 610 J = 1,3
                        STRESS(J,I) = STRESS(J,I) +
     .                                DX(J) * DFA(I) / VOLUME
  610                 CONTINUE
                    ENDIF
  620             CONTINUE
                ENDIF
  630         CONTINUE
            ENDIF
  640     CONTINUE
  650   CONTINUE

C       Find contribution of partial-core-correction
        IF (NPCC .EQ. 1) THEN
          DO 720 IA = 1,NA
            IUA = INDXUA(IA)
            IS = ISA(IA)
            RA = RCORE( IS )
            IF (RA .GT. 0.D0) THEN
C             Loop over mesh points and subpoints inside RMAX
              DO 710 IOP = 1,MOP
                IP0 = INDEXP( IPA(IA) + IDOP(IOP) )
                IF (IP0 .GT. 0) THEN
                  DO 700 ISP = 1,NSP
                    DO 660 I = 1,3
                      DX(I) = XDOP(I,IOP) + XDSP(I,ISP) - DXA(I,IA)
  660               CONTINUE
                    R = SQRT( DOT( DX, DX, 3 ) )
                    IF (R .LT. RA) THEN
                      IP = ISP + NSP * (IP0 - 1)
                      CALL CHCORE( IS, DX, RHOP, GRRHO )
                      DO 690 ISPIN = 1,NSD
                        VXC = VSCF(IP,ISPIN) - VAUX(IP)
                        DO 680 I = 1,3
                          DFA(I) = DVOL * VXC * GRRHO(I) / NSD
                          IF (IFA .EQ. 1) FA(I,IUA) = FA(I,IUA) + DFA(I)
                          IF (ISTR .EQ. 1) THEN
                            DO 670 J = 1,3
                              STRESS(J,I) = STRESS(J,I) +
     .                                      DX(J) * DFA(I) / VOLUME
  670                       CONTINUE
                          ENDIF
  680                   CONTINUE
  690                 CONTINUE
                    ENDIF
  700             CONTINUE
                ENDIF
  710         CONTINUE
            ENDIF
  720     CONTINUE
        ENDIF

C       VAUX is (minus) the potential which multiplies RHOATM
        IF (NSD .EQ. 2) THEN
          DO 730 IP = 1,NTP
            VAUX(IP) = 0.5D0 * VAUX(IP)
  730     CONTINUE
        ENDIF

C       Initialize auxiliary array for DFSCF
        DO 740 I = 1, NORB
          AUX(I) = 0.D0
  740   CONTINUE

C       Loop over orbitals
        DO 910 IO = 1,NORB

C         Find gradient of orbital IO at mesh points
          IA = IAORB(IO)
          IUA = INDXUA(IA)
          IPHI = IPHORB(IO)
          IS = ISA(IA)
          R2O = RCUT(IS,IPHI)**2
C         Loop over mesh points inside RMAX
          JOP = 0
          DO 810 IOP = 1,MOP
            IP0 = INDEXP( IPA(IA) + IDOP(IOP) )
            IF (IP0 .GT. 0) THEN
C             Loop over sub-points to find if point is within range
              WITHIN = .FALSE.
              DO 755 ISP = 1,NSP
                DO 750 I = 1,3
                  DXSP(I,ISP) = XDOP(I,IOP) + XDSP(I,ISP) - DXA(I,IA)
  750           CONTINUE
                R2SP(ISP) = DXSP(1,ISP)**2 + DXSP(2,ISP)**2 +
     .                      DXSP(3,ISP)**2
                IF (R2SP(ISP) .LT. R2O) WITHIN = .TRUE.
  755         CONTINUE
              IF (WITHIN) THEN
                JOP = JOP + 1
                N = ENDPH(IO-1) + JOP
                LISTPH(N) = IP0
                DO 800 ISP = 1,NSP
                  IF (R2SP(ISP) .LT. R2O) THEN
                    CALL PHIATM( IS, IPHI, DXSP(1,ISP),
     .                           PHIP, GRPHI(1,ISP,JOP) )
*                   PHI(ISP,N) = PHIP
                    IF (ISTR .EQ. 1) THEN
                      DO 770 I = 1,3
                        DO 760 J = 1,3
                          XGRPHI(J,I,ISP,JOP) = DXSP(J,ISP) * 
     .                                     GRPHI(I,ISP,JOP) / VOLUME
  760                   CONTINUE
  770                 CONTINUE
                    ENDIF
                  ELSE
*                   PHI(ISP,N) = 0.D0
                    DO 790 I = 1,3
                      GRPHI(I,ISP,JOP) = 0.D0
                      DO 780 J = 1,3
                        XGRPHI(J,I,ISP,JOP) = 0.D0
  780                 CONTINUE
  790               CONTINUE
                  ENDIF
  800           CONTINUE
              ENDIF
            ENDIF
  810     CONTINUE

C         Find contribution of orbital IO tO forces and stress
          DO 900 ISPIN = 1,NSPIN
            IF (IFA .EQ. 1) THEN
              CALL DFSCF( NORB, INDXUO, ENDPH,  LISTPH, PHI, IO,
     .                    3, GRPHI, NSP, NSP,
     .                    NMP,  ENDPHT, LSTPHT, INDPHT,
     .                    MAXND, NUMD, LISTD, DSCF(1,1,ISPIN), DATM,
     .                    VOLUME, VSCF(1,ISPIN), VAUX,
     .                    DFASCF, DFAATM, MAXAUX, AUX )
              DO 850 I = 1,3
                IF (ISPIN .LE. 2) THEN
                  FA(I,IUA) = FA(I,IUA) + DFASCF(I) - DFAATM(I)
                ELSE
C                 Factor 2 takes into account that there are two 
C                 nondiagonal elements of the (non-colinear) spin  
C                 density matrix, whose real and imaginary parts
C                 are stored in ISPIN=3,4
                  FA(I,IUA) = FA(I,IUA) + 2.D0 * DFASCF(I)
                ENDIF
  850         CONTINUE
            ENDIF
            IF (ISTR .EQ. 1) THEN
              CALL DFSCF( NORB, INDXUO, ENDPH,  LISTPH, PHI, IO,
     .                    9, XGRPHI, NSP, NSP,
     .                    NMP,  ENDPHT, LSTPHT, INDPHT,
     .                    MAXND, NUMD, LISTD, DSCF(1,1,ISPIN), DATM,
     .                    VOLUME, VSCF(1,ISPIN), VAUX,
     .                    DSTRSC, DSTRAT, MAXAUX, AUX )
              DO 890 J = 1,3
                DO 880 I = 1,3
                  IF (ISPIN .LE. 2) THEN
                    STRESS(I,J)= STRESS(I,J)+DSTRSC(I,J)-DSTRAT(I,J)
                  ELSE
                    STRESS(I,J)= STRESS(I,J) + 2.D0*DSTRSC(I,J)
                  ENDIF
  880           CONTINUE
  890         CONTINUE
            ENDIF
  900     CONTINUE
  910   CONTINUE

C       Stop time counter for force calculation part
        CALL TIMER( 'DHSCF4', 2 )
      ENDIF
C ---------------------------------------------------------------------

C Stop time counter -----------------------------------------------
  999 CONTINUE
      CALL TIMER( 'DHSCF', 2 )
C -----------------------------------------------------------------

      FRSTME = .FALSE.
      END




