
      PROGRAM PLOTDEN
C **********************************************************************
C Reads density matrix from siesta and calculates the density of charge 
C at the points of a plane in real space 
C Coded by J. Junquera 11/98
C
C Version: 0.1.1
C **********************************************************************

      IMPLICIT NONE

      INCLUDE 'siesta.h'
      INCLUDE 'fdfdefs.h'

      INTEGER 
     .  NUO, NO, NA, NSPIN, 
     .  LASTO(0:MAXA), ISA(MAXA), INDXUO(MAXO),
     .  IPHORB(MAXO),
     .  NUMH(MAXO), LISTH(MAXNO,MAXO)

      INTEGER
     .  IOPTION, NPX, NPY, IUNITCD 

      DOUBLE PRECISION
     .  CELL(3,3), XA(3,MAXA), RMAXO, VOLUME, VOLCEL

      DOUBLE PRECISION
     .  DSCF(MAXNO,MAXO,MAXSPN), DSCFNA(MAXNO,MAXO,MAXSPN),
     .  DATM(MAXO)

      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX, COORPO(3,3), NORMAL(3), DIRVER1(3),
     .  DIRVER2(3), ARMUNI

      LOGICAL 
     .  FOUND

      CHARACTER
     .  FILEIN*20, FILEOUT*20

      EXTERNAL
     .  FDF_INIT, REDATA, VOLCEL, READPLA, IODM, DMNA, RHOOFR 

      DATA NORMAL /0.D0,0.D0,1.D0/
      DATA COORPO /1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0/
      DATA DIRVER1 /1.D0,0.D0,0.D0/
      DATA DIRVER2 /0.D0,1.D0,0.D0/


C **********************************************************************
C INTEGER NUO                      : Number of orbitals in unit cell
C INTEGER NO                       : Number of atomic orbitals in supercell
C INTEGER NA                       : Number of atoms in supercell
C INTEGER NSPIN                    : Number of different spin polarizations
C                                    Nspin = 1 => Unpolarized, Nspin = 2 => Pol.
C INTEGER LASTO(0:MAXA)            : Position of last orbital of each atom
C INTEGER ISA(MAXA)                : Species index of each atom
C INTEGER INDXUO(MAXO)             : Equivalent orbital in unit cell
C INTEGER IPHORB(MAXO)             : Orbital index of each orbital in its atom
C INTEGER NUMH(MAXO)               : Control vector of Density Matrix
C                                    (number of nonzero elements of each row)
C INTEGER LISTH(MAXNO,MAXO)        : Control vector of Density Matrix
C                                    (list of nonzero elements of each row)
C INTEGER IUNITCD                  : Units of the charge density
C REAL*8  CELL(3,3)                : Unit cell vectors CELL(IXYZ,IVECT)
C REAL*8  VOLUME                   : Volumen of unit cell (in bohr**3)
C REAL*8  XA(3,MAXA)               : Atomic positions in cartesian coordinates
C REAL*8  RMAXO                    : Maximum range of basis orbitals
C REAL*8  DSCF(MAXNO,MAXO,MAXSPN)  : Density Matrix (DM)
C REAL*8  DSCFNA(MAXNO,MAXO,MAXSPN): Density Matrix for Neutral Atoms
C REAL*8  DATM(MAXO)               : Neutral atom charge of each orbital
C LOGICAL FOUND                    : Has DM been found in disk?
C                                    (Only when task = 'read')
C INTEGER IOPTION                  : Option to generate the plane
C                                    1 = Normal vector
C                                    2 = Two vectors belonging to the plane
C                                    3 = Three points of the plane
C INTEGER NPX, NPY                 : Number of points generated along x and y
C                                    directions ina a system of reference 
C                                    in which the third component of the 
C                                    points of the plane is zero 
C                                    (Plane Reference Frame; PRF)
C REAL*8  XMIN, XMAX               : Limits of the plane in the PRF 
C                                    for x-direction
C REAL*8  YMIN, YMAX               : Limits of the plane in the PRF 
C                                    for y-direction
C REAL*8  NORMAL(3)                : Components of the normal vector 
C                                    used to define the plane
C REAL*8  DIRVER(3)                : Components of the two vectors contained 
C                                    in the plane
C REAL*8  COORPO(3,3)              : Coordinates of the three points used 
C                                    to define the plane
C **********************************************************************

C Set up fdf -----------------------------------------------------------
      FILEIN  = 'stdin'
      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

C Read tables calculated in siesta from files --------------------------
      CALL REDATA( NSPIN, NUO, NO, NA, 
     .             MAXO, MAXA, MAXNO, MAXSPN, 
     .             CELL, RMAXO, XA,
     .             LASTO, ISA, IPHORB, DATM,
     .             NUMH, LISTH, INDXUO )

C Calculate the volumen of the unit cell -------------------------------
      VOLUME = VOLCEL( CELL )

C Read Density Matrix from files ---------------------------------------
      CALL IODM('READ', MAXNO, MAXO, NUO, NSPIN,
     .          NUMH, LISTH, DSCF, FOUND )
      IF (.NOT. FOUND) THEN
        WRITE(6,*)' DENSITY MATRIX NOT FOUND              '
        WRITE(6,*)' CHECK YOU HAVE COPY IT FROM THE       '
        WRITE(6,*)' DIRECTORY WHERE YOU HAVE RUN SIESTA   '
        STOP
      ENDIF 

C Read option to generate the plane ------------------------------------
      CALL READPLA( IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .              NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2,
     .              ARMUNI, MAXA, XA, VOLUME, IUNITCD )

cC Dump input inot the output files -------------------------------------
c      CALL WROUT( NSPIN, IOPTION)

C Form Density Matrix for Neutral and Isolated Atoms -------------------
      CALL DMNA( NO, NSPIN, MAXO, MAXNO, NUMH, LISTH, INDXUO, 
     .           DATM, DSCFNA )

C Calulate the charge density ------------------------------------------
      CALL RHOOFR( CELL, NA, MAXNA, MAXNO, MAXA, MAXO, MAXSPN,
     .             NSPIN, RMAXO, LASTO, XA, ISA, IPHORB,
     .             INDXUO, NUMH, LISTH, DSCF, DSCFNA,
     .             IOPTION, XMIN, XMAX, YMIN, YMAX,
     .             NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2,
     .             ARMUNI, IUNITCD )

      END
