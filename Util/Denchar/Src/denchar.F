
      PROGRAM DENCHAR

C **********************************************************************
C Reads density matrix from SIESTA and calculates the charge density
C at the points of a plane in real space 
C Coded by J. Junquera 11/98
C Modified by J. Junquera 07/01 
C
C Version: 1.0.1
C **********************************************************************
C
C  Modules
C
      USE PRECISION
      USE PARALLEL
      USE BASIS_IO
      USE LISTSC_MODULE, ONLY: LISTSC_INIT
      USE FDF
#ifdef MPI
      use mpi
#endif

      IMPLICIT NONE

      INTEGER
     .   NO_U, NO_L, NO_S, NA_S, NSPIN, MAXND, MAXNA,
     .   NODE, NODES, NSC(3)

      INTEGER
     .  IOPTION, NPX, NPY, IUNITCD, ISCALE 

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  ISA, LASTO, IPHORB, INDXUO, 
     .  NUMD, LISTD, LISTDPTR

#ifdef MPI
      INTEGER
     .  MPIerror 
#endif

      DOUBLE PRECISION
     .   CELL(3,3), VOLUME, VOLCEL, RMAXO

      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX, COORPO(3,3), NORMAL(3), DIRVER1(3),
     .  DIRVER2(3), ARMUNI

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     .   XA, DSCF

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::
     .   DATM, DSCFNA

      CHARACTER
     .  FILEIN*20, FILEOUT*20

      LOGICAL 
     .  FOUND

      EXTERNAL
     .  DMNA, IODM, READPLA, REDATA, REINIT, RHOOFR, VOLCEL

      DATA NORMAL /0.D0,0.D0,1.D0/
      DATA COORPO /1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0/
      DATA DIRVER1 /1.D0,0.D0,0.D0/
      DATA DIRVER2 /0.D0,1.D0,0.D0/

C ****** READ FROM SIESTA **********************************************
C INTEGER NO_U                : Total number of orbitals in the unit cell
C INTEGER NO_L                : Number of orbitals stored in a local node
C INTEGER NO_S                : Total number of orbitals in the supercell
C INTEGER NA_S                : Total number of atoms in the supercell
C INTEGER NSPIN               : Number of different spin polarizations
C                               Nspin = 1 => Unpolarized, Nspin = 2 => Polarized
C INTEGER MAXND               : Maximum number
C                               of basis orbitals interacting, either directly
C                               or through a KB projector, with any orbital
C INTEGER MAXNA               : Maximum number of neighbours of any atom
C INTEGER NODE                : The local processor
C INTEGER NODES               : The total number of processors
C INTEGER NSC(3)              : Num. of unit cells in each supercell direction
C INTEGER ISA(MAXA)           : Species index of each atom in the supercell
C INTEGER LASTO(0:MAXA)       : Last orbital of each atom in array iphorb
C INTEGER IPHORB(MAXO)        : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)        : Equivalent orbital in unit cell
C INTEGER NUMD(NO_U)          : Number of nonzero elements of each row of the
C                               Hamiltonian matrix between atomic orbitals
C INTEGER LISTD(MAXND)        : Nonzero Hamiltonian-matrix element
C                               column indexes for each matrix row
C                               For parallel execution, listh contains the
C                               elements for rows that involve
C                               any locally stored orbitals. In the case
C                               where parallelisation is over K points then
C                               the full listh matrix is needed on every
C                               Node.
C INTEGER LISTDPTR(NO_U)      : Pointer to where each row of listh starts - 1
C                               The reason for pointing to the element before
C                               the first one is so that when looping over the
C                               elements of a row there is no need to shift by
C                               minus one.
C REAL*8  CELL(3,3)           : Supercell vectors CELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  VOLUME              : Volumen of unit cell (in bohr**3)
C REAL*8  RMAXO               : Maximum range of basis orbitals
C REAL*8  XA(3,NA_S)          : Atomic coordinates in cartesian coordinates
C                               (units in bohrs)
C REAL*8  DATM(NO_S)          : Occupations of basis orbitals in free atom
C REAL*8  DSCF(MAXND,NSPIN)   : Density Matrix (DM)
C ****** INFORMATION OF THE PLANE **************************************
C INTEGER IOPTION             : Option to generate the plane
C                               1 = Normal vector
C                               2 = Two vectors belonging to the plane
C                               3 = Three points of the plane
C                               4 = Three atomic indices
C INTEGER NPX, NPY            : Number of points generated along x and y
C                               directions in a system of reference 
C                               in which the third component of the 
C                               points of the plane is zero 
C                               (Plane Reference Frame; PRF)
C INTEGER IUNITCD             : Units of the charge density
C INTEGER ISCALE              : Units of the points of the plane
C REAL*8  XMIN, XMAX          : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX          : Limits of the plane in the PRF for y-direction
C REAL*8  COORPO(3,3)         : Coordinates of the three points used 
C                               to define the plane
C REAL*8  NORMAL(3)           : Components of the normal vector 
C                               used to define the plane
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


C Set up fdf -----------------------------------------------------------
      FILEIN  = 'stdin'
      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

C Read some variables from SIESTA to define the limits of some arrays --
      CALL REINIT( NO_S, NA_S, NO_U, MAXND, MAXNA, NSPIN )

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
      CALL REDATA( NO_S, NA_S, NO_U, MAXND, NSPIN,
     .             ISA, IPHORB, INDXUO, LASTO,
     .             CELL, NSC, XA, RMAXO, DATM )

C Read the information about the basis set -----------------------------
      CALL READ_BASIS_ASCII

C Initialize listsc ----------------------------------------------------
      CALL LISTSC_INIT( NSC, NO_U )

C Calculate the volume of the unit cell --------------------------------
      VOLUME = VOLCEL( CELL )

C Initialise MPI and set processor number ------------------------------
#ifdef MPI
      CALL MPI_Init( MPIerror )
      CALL MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      CALL MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      NODE =  0 
      NODES = 1
#endif

C Find number of locally stored orbitals and allocated related arrays --
      CALL GetNodeOrbs(NO_U,NODE,NODES,NO_L)

      ALLOCATE(LISTDPTR(NO_L))
      CALL MEMORY('A','I',NO_L,'denchar')
      LISTDPTR(:) = 0

      ALLOCATE(NUMD(NO_L))
      CALL MEMORY('A','I',NO_L,'denchar')
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
      CALL IODM('READ', MAXND, NO_L, NSPIN,
     .          NUMD, LISTDPTR, LISTD, DSCF, FOUND )
      IF (.NOT. FOUND) THEN
        WRITE(6,*)' DENSITY MATRIX NOT FOUND              '
        WRITE(6,*)' CHECK YOU HAVE COPY IT FROM THE       '
        WRITE(6,*)' DIRECTORY WHERE YOU HAVE RUN SIESTA   '
        STOP
      ENDIF 

C Read option to generate the plane ------------------------------------
      CALL READPLA( NA_S, XA, VOLUME, 
     .              IOPTION, IUNITCD, ISCALE, NPX, NPY, 
     .              XMIN, XMAX, YMIN, YMAX, 
     .              COORPO, NORMAL, DIRVER1, DIRVER2,
     .              ARMUNI )

C Form Density Matrix for Neutral and Isolated Atoms -------------------
      CALL DMNA( NO_U, NO_S, MAXND, NUMD, LISTD, LISTDPTR,
     .           DATM, DSCFNA )
      
      CALL RHOOFR( NA_S, NO_S, NO_L, MAXND, MAXNA, NSPIN, 
     .             ISA, IPHORB, INDXUO, LASTO,
     .             XA, CELL, NUMD, LISTD, LISTDPTR, DSCF, DSCFNA,
     .             IOPTION, XMIN, XMAX, YMIN, YMAX,
     .             NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2,
     .             ARMUNI, IUNITCD, ISCALE, RMAXO )


      END PROGRAM DENCHAR
