! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      PROGRAM MAINSTM

C **********************************************************************
C Reads eigenfncts and eigenvalues from SIESTA and computes STM images
C
C Coded by P. Ordejon, November 2004
C
C Modified by N. Lorente, August 2005
C
C **********************************************************************
C
C  Modules
C
      USE PRECISION
      USE BASIS_IO
      USE LISTSC_MODULE, ONLY: LISTSC_INIT
      USE FDF
      use m_wfsx

      IMPLICIT NONE

      INTEGER ns,
     .   NO_U, NO_S, NA_S, NSPIN, MAXND, MAXNA,
     .   NSC(3), NK, NUMWF

      INTEGER
     .  NPX, NPY, NPZ, IUNITCD, IX, IV, NSCX, NSCY

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  ISA, LASTO, IPHORB, INDXUO, 
     .  NUMD, LISTD, LISTDPTR, NWF

      real(dp), DIMENSION(:), ALLOCATABLE  ::  WK
      INTEGER, DIMENSION(:,:), ALLOCATABLE ::  INDW

      REAL(DP)
     .   CELL(3,3), VOLUME, VOLCEL, RMAXO, UCELL(3,3), V0, EMAX, EMIN

      REAL(DP)  ZMIN, ZMAX, ZREF, ARMUNI

      REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE ::  RPSI,IPSI
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE ::   E
      REAL(DP), DIMENSION(:,:), ALLOCATABLE ::   XA, K

      CHARACTER  FILEIN*20, FILEOUT*20, SNAME*40, FNAME*50

      INTEGER :: wf_unit
      logical :: gamma, non_coll
      integer :: nspin_blocks
      
      EXTERNAL :: IODM, REDATA, REINITSTM, VOLCEL


C ****** READ FROM SIESTA **********************************************
! integer ns                  : Total number of species
C INTEGER NO_U                : Total number of orbitals in the unit cell
C INTEGER NO_S                : Total number of orbitals in the supercell
C INTEGER NA_S                : Total number of atoms in the supercell
C INTEGER NSPIN               : Number of different spin polarizations
C                               Nspin = 1 => Unpolarized, Nspin = 2 => Polarized
C INTEGER MAXND               : Maximum number
C                               of basis orbitals interacting, either directly
C                               or through a KB projector, with any orbital
C INTEGER MAXNA               : Maximum number of neighbours of any atom
C INTEGER NSC(3)              : Num. of unit cells in each supercell direction
C INTEGER NUMWF               : Max num of wavefncts to print for a given k-po.
C INTEGER NWF(NK)             : Num of wavefncts to print for each k-point
C REAL(dp) WK(NK)             : weight of each k-point
C INTEGER ISA(MAXA)           : Species index of each atom in the supercell
C INTEGER LASTO(0:MAXA)       : Last orbital of each atom in array iphorb
C INTEGER IPHORB(MAXO)        : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)        : Equivalent orbital in unit cell
C INTEGER NUMD(NO_U)          : Number of nonzero elements of each row of the
C                               Hamiltonian matrix between atomic orbitals
C INTEGER LISTD(MAXND)        : Nonzero Hamiltonian-matrix element
C                               column indexes for each matrix row
C INTEGER LISTDPTR(NO_U)      : Pointer to where each row of listh starts - 1
C                               The reason for pointing to the element before
C                               the first one is so that when looping over the
C                               elements of a row there is no need to shift by
C                               minus one.
C INTEGER NK                  : Number of k-points in wave functions file
C INTEGER INDW(NK,NUMWF)      : Indices of wavefuncs to print for each k-point
C REAL*8  CELL(3,3)           : Supercell vectors CELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  UCELL(3,3)          : Unit cell vectors UCELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  VOLUME              : Volumen of unit cell (in bohr**3)
C REAL*8  RMAXO               : Maximum range of basis orbitals
C REAL*8  XA(3,NA_S)          : Atomic coordinates in cartesian coordinates
C                               (units in bohrs)
C REAL*8  RPSI(NO_U,NK,NUMWF,NSPIN): Wave function coefficients (real part)
C REAL*8  IPSI(NO_U,NK,NUMWF,NSPIN): Wave function coefficients (imag part)
C REAL*8  E(NK,NUMWF,NSPIN)        : Wave function energies
C REAL*8  K(NK,3)             : K-points
C ****** INFORMATION OF THE PLANE OR 3D GRID ***************************
C INTEGER NPX, NPY, NPZ       : Number of points generated along x and y
C                               (and z for 3D grids) directions in a system of 
C                               reference in which the third component of the 
C                               points of the plane is zero 
C                               (Plane Reference Frame; PRF)
C INTEGER IUNITCD             : Units of the charge density
C REAL*8  ZREF                : Position of reference plane for wf. extrapol.
C REAL*8  ZMIN, ZMAX          : Limits of the or 3D grid in the PRF for z-axis
C INTEGER NSCX, NSCY          : Number of cells to plot in cube file (in x
C                               and y directions)
C REAL*8  ARMUNI              : Conversion factors for the charge density
C **********************************************************************


C Set up fdf -----------------------------------------------------------
      FILEIN  = 'stdin'
      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

      WRITE(6,*)
      WRITE(6,*) 'STM Simulation program '
      WRITE(6,*) 'P. Ordejon and N. Lorente, Nov. 04'
      WRITE(6,*)
      WRITE(6,*) 'Reading information from SIESTA'
      WRITE(6,*)


C Read some variables from SIESTA to define the limits of some arrays --
      CALL REINITSTM( NO_S, NA_S, NO_U, MAXND, MAXNA, NSPIN)

C Allocate some variables ----------------------------------------------
      ALLOCATE(XA(3,NA_S))
      CALL MEMORY('A','D',3*NA_S,'mainstm')

      ALLOCATE(LASTO(0:NA_S))
      CALL MEMORY('A','D',NA_S+1,'mainstm')

      ALLOCATE(ISA(NA_S))
      CALL MEMORY('A','D',NA_S,'mainstm')

      ALLOCATE(IPHORB(NO_S))
      CALL MEMORY('A','D',NO_S,'mainstm')

      ALLOCATE(INDXUO(NO_S))
      CALL MEMORY('A','D',NO_S,'mainstm')


C Read some variables from SIESTA --------------------------------------
      CALL REDATA( NO_S, NA_S, NO_U, MAXND, NSPIN,
     .             ISA, IPHORB, INDXUO, LASTO,
     .             CELL, NSC, XA, RMAXO )

C Read the information about the basis set -----------------------------
      CALL READ_BASIS_ASCII(ns)

C Initialize listsc ----------------------------------------------------
      CALL LISTSC_INIT( NSC, NO_U )

C Calculate the volume of the unit cell --------------------------------
      VOLUME = VOLCEL( CELL ) / (NSC(1) * NSC(2) * NSC(3))

C Calculate unit cell vectors
      DO IX=1,3
        DO IV=1,3
          UCELL(IX,IV) = CELL(IX,IV)/NSC(IV)
        ENDDO
      ENDDO


C Allocate variables

C call readwaves just to find out number of wavefunctions to print
C allocate temporary space for eigenvalues and eigenfunctions
      
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = trim(SNAME) //'.WFSX'

      CALL IO_ASSIGN(wf_unit)
      OPEN (wf_unit, FILE=FNAME, FORM='unformatted',
     $       STATUS='unknown',position='rewind')

      CALL wfsx_get_info(wf_unit,NSPIN,NO_U,NUMWF,NK,gamma)

      non_coll = (nspin == 4)
      nspin_blocks = nspin
      if (non_coll) nspin_blocks = 1

      print *, "NK, NUMWF, gamma: ", NK, NUMWF, gamma
      
      ! Prepare file to stream over it in routine stm
      
      rewind(wf_unit)
      read(wf_unit)
      read(wf_unit)
      read(wf_unit)
      read(wf_unit)             ! advance past header section
        

C Read option for STM 
      CALL READSTM( VOLUME, 
     .     IUNITCD, NPX, NPY, NPZ, ZREF, ZMIN, ZMAX, EMAX, EMIN,
     .      NSCX, NSCY, ARMUNI )

C Call routine to calculate the value of the potential at vacuum

      CALL VACPOT(V0)

C Call routine to compute STM image

        CALL STM   ( NA_S, NO_S, NO_U, MAXNA, nspin_blocks, non_coll,
     .               ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .               wf_unit, NK, gamma,
     .               ZREF, ZMIN, ZMAX, NPX, NPY, NPZ, NSCX, NSCY,
     .               V0, EMAX, EMIN,
     .               ARMUNI, IUNITCD, RMAXO )


      END PROGRAM MAINSTM
