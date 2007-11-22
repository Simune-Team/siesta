! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      PROGRAM DENCHAR

C **********************************************************************
C Reads density matrix from SIESTA and calculates the charge density
C at the points of a plane in real space, or at a 3D grid of points
C Coded by J. Junquera 11/98
C Modified by J. Junquera 07/01 
C 3D and wavefunction capabilities coded by P. Ordejon, June 2003
C
C Version: 1.2.0
C **********************************************************************
C
C  Modules
C
      use f2kcli
      use m_denchar_init
      use m_denchar_io
      use m_denchar_work
      use m_denchar_geom

      USE PRECISION
      USE PARALLEL
      USE BASIS_IO
      USE LISTSC_MODULE, ONLY: LISTSC_INIT
      USE FDF
      use sys, only: die

C     It has to be MPI-aware since it calls routines down the line...
C
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

#ifdef MPI
      integer  MPIerror
#endif

      INTEGER
     .   NO_U, NO_S, NA_S, NSPIN, MAXND, MAXNA,
     .   NSC(3), NWF

      INTEGER
     .  IDIMEN, IOPTION, NPX, NPY, NPZ, IUNITCD, ISCALE 

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  ISA, LASTO, IPHORB, INDXUO, 
     .  NUMD, LISTD, LISTDPTR, INDW


      DOUBLE PRECISION
     .   CELL(3,3), VOLUME, VOLCEL, RMAXO

      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3),
     .  DIRVER2(3), ARMUNI

      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::
     .   PSI

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     .   XA, DSCF, E

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::
     .   DATM, DSCFNA

      CHARACTER
     .  FILEIN*20, FILEOUT*20

      LOGICAL 
     .  FOUND, CHARGE, WAVES

      EXTERNAL  IODM, VOLCEL

      integer narg, ns

      DATA NORMAL /0.D0,0.D0,1.D0/
      DATA COORPO /1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0/
      DATA DIRVER1 /1.D0,0.D0,0.D0/
      DATA DIRVER2 /0.D0,1.D0,0.D0/

C ****** READ FROM SIESTA **********************************************
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
C INTEGER NWF                 : Number of wavefunctions to print
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
C INTEGER INDW(NWF)           : Indices of wavefunctions to print
C REAL*8  CELL(3,3)           : Supercell vectors CELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  VOLUME              : Volumen of unit cell (in bohr**3)
C REAL*8  RMAXO               : Maximum range of basis orbitals
C REAL*8  XA(3,NA_S)          : Atomic coordinates in cartesian coordinates
C                               (units in bohrs)
C REAL*8  DATM(NO_S)          : Occupations of basis orbitals in free atom
C REAL*8  DSCF(MAXND,NSPIN)   : Density Matrix (DM)
C REAL*8  PSI(NO_U,NWF,NSPIN) : Wave function coefficients
C REAL*8  E(NWF,NSPIN)        : Wave function energies
C ****** INFORMATION OF THE PLANE OR 3D GRID ***************************
C INTEGER IDIMEN              : Specifies 2D or 3D mode
C LOGICAL CHARGE              : Should charge density be computed?
C LOGICAL WAVES               : Should wave functions be computed?
C INTEGER IOPTION             : Option to generate the plane or 3D grid
C                               1 = Normal vector to xy plane (ie, z direction)
C                               2 = Two vectors belonging to the xy plane
C                               3 = Three points of the xy plane
C                               4 = Three atomic indices define the xy plane
C INTEGER NPX, NPY, NPZ       : Number of points generated along x and y
C                               (and z for 3D grids) directions in a system of 
C                               reference in which the third component of the 
C                               points of the plane is zero 
C                               (Plane Reference Frame; PRF)
C INTEGER IUNITCD             : Units of the charge density
C INTEGER ISCALE              : Units of the points of the plane or 3D grid
C REAL*8  XMIN, XMAX          : Limits of the plane in the PRF for x-axis
C REAL*8  YMIN, YMAX          : Limits of the plane in the PRF for y-axis
C REAL*8  ZMIN, ZMAX          : Limits of the or 3D grid in the PRF for z-axis
C REAL*8  COORPO(3,3)         : Coordinates of the three points used 
C                               to define the xy plane
C REAL*8  NORMAL(3)           : Components of the normal vector 
C                               used to define the xy plane (z-direction)
C REAL*8  DIRVER1(3)          : Components of the first vector contained 
C                               in the plane
C REAL*8  DIRVER2(3)          : Components of the first vector contained 
C                               in the plane
C REAL*8  ARMUNI              : Conversion factors for the charge density
C ****** INTERNAL VARIABLES ********************************************
C REAL*8  DSCFNA(MAXND)       : Density Matrix for Neutral Atoms
C LOGICAL FOUND               : Has DM been found in disk?
C                               (Only when task = 'read')
C **********************************************************************


C Initialise MPI and set processor number
#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node =  0
      Nodes = 1
#endif
      if (Nodes /= 1)
     $            print *, "You are running denchar on several nodes!"

      IONode = (Node .eq. 0)

C Set up fdf -----------------------------------------------------------
!      narg = command_argument_count()
!      if (narg == 0) then
!         FILEIN = "stdin"
!      else if (narg == 1) then
!         call get_command_argument(1,FILEIN)
!      else
!         stop "too many command-line arguments"
!      endif
      FILEIN='stdin'

      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

C Read some variables from SIESTA to define the limits of some arrays --
      CALL REINIT( NO_S, NA_S, NO_U, MAXND, MAXNA, NSPIN, IDIMEN,
     .            CHARGE, WAVES )

C Allocate some variables ----------------------------------------------
      ALLOCATE(XA(3,NA_S))
      CALL MEMORY('A','D',3*NA_S,'denchar')

      ALLOCATE(LASTO(0:NA_S))
      CALL MEMORY('A','D',NA_S+1,'denchar')

      ALLOCATE(ISA(NA_S))
      CALL MEMORY('A','D',NA_S,'denchar')

      ALLOCATE(IPHORB(NO_S))
      CALL MEMORY('A','D',NO_S,'denchar')

      ALLOCATE(INDXUO(NO_S))
      CALL MEMORY('A','D',NO_S,'denchar')

      ALLOCATE(DATM(NO_S))
      CALL MEMORY('A','D',NO_S,'denchar')

C Read some variables from SIESTA --------------------------------------
      CALL REDATA( NO_S, NA_S, 
     .             ISA, IPHORB, INDXUO, LASTO,
     .             CELL, NSC, XA, RMAXO, DATM )

C Read the information about the basis set -----------------------------
      CALL READ_BASIS_ASCII(ns)

C Initialize listsc ----------------------------------------------------
      CALL LISTSC_INIT( NSC, NO_U )

C Calculate the volume of the unit cell --------------------------------
      VOLUME = VOLCEL( CELL ) / (NSC(1) * NSC(2) * NSC(3))


C Allocate variables

C If this is a charge calculation, allocate space for DM
      IF (CHARGE) THEN
        ALLOCATE(LISTDPTR(NO_U))
        CALL MEMORY('A','I',NO_U,'denchar')
        LISTDPTR(:) = 0

        ALLOCATE(NUMD(NO_U))
        CALL MEMORY('A','I',NO_U,'denchar')
        NUMD(:) = 0

C Allocate some other variables ----------------------------------------
        IF (.NOT.ALLOCATED(LISTD)) THEN
          ALLOCATE(LISTD(MAXND))
          CALL MEMORY('A','I',MAXND,'denchar')
        ENDIF
        
        IF (ALLOCATED(DSCF)) THEN
          CALL MEMORY('D','D',SIZE(DSCF),'denchar')
          DEALLOCATE(DSCF)
        ENDIF
        ALLOCATE(DSCF(MAXND,NSPIN))
        CALL MEMORY('A','D',MAXND*NSPIN,'denchar')

        IF (ALLOCATED(DSCFNA)) THEN
          CALL MEMORY('D','D',SIZE(DSCFNA),'denchar')
          DEALLOCATE(DSCFNA)
        ENDIF
        ALLOCATE(DSCFNA(MAXND))
        CALL MEMORY('A','D',MAXND,'denchar')
  
C Read Density Matrix from files ---------------------------------------
        CALL IODM('READ', MAXND, NO_U, NSPIN,
     .            NUMD, LISTDPTR, LISTD, DSCF, FOUND )
        IF (.NOT. FOUND) THEN
          WRITE(6,*)' DENSITY MATRIX NOT FOUND              '
          WRITE(6,*)' CHECK YOU HAVE COPIED IT FROM THE       '
          WRITE(6,*)' DIRECTORY WHERE YOU HAVE RUN SIESTA   '
          call die()
        ENDIF 
      ENDIF

      IF (WAVES) THEN
C call readwaves just to find out number of wavefunctions to print
C allocate temporary space for eigenvalues and eigenfunctions
        NWF = 1
        ALLOCATE(INDW(NWF))
        CALL MEMORY('A','I',NWF,'denchar')
        ALLOCATE(E(NWF,NSPIN))
        CALL MEMORY('A','D',NWF*NSPIN,'denchar')
        ALLOCATE(PSI(NO_U,NWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NWF*NSPIN,'denchar')

        CALL READWAVES(NSPIN,NO_U,0,NWF,PSI,E,INDW)

C deallocate temporary space 
        CALL MEMORY('D','I',SIZE(INDW),'denchar')
        DEALLOCATE(INDW)
        CALL MEMORY('D','D',SIZE(E),'denchar')
        DEALLOCATE(E)
        CALL MEMORY('D','D',SIZE(PSI),'denchar')
        DEALLOCATE(PSI)

C allocate space for eigenenergies and eigenfunctions, assumed to be real
        ALLOCATE(INDW(NWF))
        CALL MEMORY('A','I',NWF,'denchar')
        ALLOCATE(E(NWF,NSPIN))
        CALL MEMORY('A','D',NWF*NSPIN,'denchar')
        ALLOCATE(PSI(NO_U,NWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NWF*NSPIN,'denchar')

C call readwaves again to actually read wavefunctions
        CALL READWAVES(NSPIN,NO_U,1,NWF,PSI,E,INDW)
      ENDIF



C Read option to generate the plane or 3D-grid -------------------------
      CALL READPLA( NA_S, XA, VOLUME, IDIMEN,
     .              IOPTION, IUNITCD, ISCALE, NPX, NPY, NPZ,
     .              XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .              COORPO, NORMAL, DIRVER1, DIRVER2, 
     .              ARMUNI )

C Form Density Matrix for Neutral and Isolated Atoms -------------------
      IF (CHARGE) THEN
        CALL DMNA( NO_U, NO_S, MAXND, NUMD, LISTD, LISTDPTR,
     .             DATM, DSCFNA )
      
        CALL RHOOFR( NA_S, NO_S, NO_U, MAXND, MAXNA, NSPIN, 
     .               ISA, IPHORB, INDXUO, LASTO,
     .               XA, CELL, NUMD, LISTD, LISTDPTR, DSCF, DSCFNA,
     .               IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .               ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .               DIRVER1, DIRVER2, 
     .               ARMUNI, IUNITCD, ISCALE, RMAXO )
      ENDIF

      IF (WAVES) THEN
        CALL WAVOFR( NA_S, NO_S, NO_U, MAXNA, NSPIN, 
     .               ISA, IPHORB, INDXUO, LASTO,
     .               XA, CELL, PSI, INDW, NWF,
     .               IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .               ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .               DIRVER1, DIRVER2, 
     .               ARMUNI, IUNITCD, ISCALE, RMAXO )
      ENDIF


      END PROGRAM DENCHAR
