
      SUBROUTINE RHOOFR( NA, NO, NUO, MAXND, MAXNA, NSPIN, 
     .                   ISA, IPHORB, INDXUO, LASTO, 
     .                   XA, CELL, NUMD, LISTD, LISTDPTR, DSCF, DSCFNA, 
     .                   IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                   NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                   ARMUNI, IUNITCD, ISCALE, RMAXO )
C **********************************************************************
C Compute the density of charge at the points of a plane in real space
C Coded by J. Junquera November'98
C **********************************************************************

      USE FDF
      USE ATMFUNCS
      USE LISTSC_MODULE, ONLY: LISTSC

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NUO, IOPTION, NPX, NPY, ISCALE, IUNITCD,
     .  MAXND, NSPIN, MAXNA,
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA),
     .  NUMD(NUO), LISTDPTR(NUO), LISTD(MAXND)

      DOUBLE PRECISION, INTENT(IN) ::
     .  XA(3,NA), XMIN, XMAX, YMIN, YMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), ARMUNI,
     .  RMAXO

      DOUBLE PRECISION, INTENT(IN) ::
     . CELL(3,3), DSCF(MAXND,NSPIN),
     . DSCFNA(MAXND)

C ****** INPUT *********************************************************
C INTEGER NA               : Total number of atoms in Supercell
C INTEGER NO               : Total number of orbitals in Supercell
C INTEGER NUO              : Total number of orbitals in Unit Cell
C INTEGER MAXND            : Maximum number
C                            of basis orbitals interacting, either directly
C                            or through a KB projector, with any orbital
C INTEGER MAXNA            : Maximum number of neighbours of any atom
C INTEGER NSPIN            : Number of different spin polarizations
C                            Nspin = 1 => unpolarized, Nspin = 2 => polarized
C INTEGER ISA(NA)          : Species index of each atom
C INTEGER IPHORB(NO)       : Orital index of each orbital in its atom
C INTEGER INDXUO(NO)       : Equivalent orbital in unit cell
C INTEGER LASTO(0:NA)      : Last orbital of each atom in array iphorb
C REAL*8  XA(3,NA)         : Atomic positions in cartesian coordinates
C                            (in bohr)
C REAL*8  CELL(3,3)        : Supercell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C INTEGER NUMD(NUO)        : Control vector of Density Matrix
C                            (number of non-zero elements of eaach row)
C INTEGER LISTDPTR(NUO)    : Pointer to where each row of listh starts - 1
C                            The reason for pointing to the element before
C                            the first one is so that when looping over the
C                            elements of a row there is no need to shift by
C                            minus one.
C INTEGER LISTD(MAXND)     : Control vector of Density Matrix
C                            (list of non-zero elements of each row)
C REAL*8  DSCF(MAXND,NSPIN): Density Matrix
C REAL*8  DSCFNA(MAXND)    : Density Matrix for Neutral Atoms
C INTEGER IOPTION          : Option to generate the plane
C                            1 = Normal vector
C                            2 = Two vectors contained in the plane
C                            3 = Three points in the plane
C INTEGER NPX,NPY          : Number of points generated along x and y 
C                            direction in a system of reference in which
C                            the third components od the points of the plane is
C                            zero (Plane Reference Frame; PRF)
C REAL*8  XMIN, XMAX       : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX       : Limits of the plane in the PRF for y-direction
C REAL*8  NORMAL(3)        : Components of the normal vector used to define 
C                            the plane
C REAL*8  DIRVER1(3)       : Components of the first vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  DIRVER2(3)       : Components of the second vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  COORPO(3,3)      : Coordinates of the three points used to define
C                            the plane (Only used if ioption = 3)
C INTEGER IUNITCD          : Unit of the charge density
C INTEGER ISCALE           : Unit if the points of the plane
C REAL*8  ARMUNI           : Conversion factor for the charge density
C REAL*8  RMAXO            : Maximum range of basis orbitals
C **********************************************************************

      INTEGER
     .  NPLAMAX, NAPLA

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  INDICES, JNA

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::
     .   R2IJ, DISCF, DIATM, DIUP, DIDOWN

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     .   PLAPO, POINRE, XAPLA, XIJ

      INTEGER
     .  NPO, IA, ISEL, NNA, UNIT1, UNIT2, UNIT3, UNIT4
 
      INTEGER
     .  I, J, IX, IN, IAT1, IAT2, JO, IO, IUO, IAVEC, IAVEC1, IAVEC2, 
     .  IS1, IS2, IPHI1, IPHI2, IND

      DOUBLE PRECISION
     .  RMAX, XPO(3), RMAX2, XVEC1(3),
     .  XVEC2(3), PHIMU, PHINU, GRPHIMU(3), GRPHINU(3)

      DOUBLE PRECISION
     .  DENCHAR, DENCHAR0, DENUP, DENDOWN
 
      LOGICAL FIRST

      CHARACTER
     .  SNAME*30, FNAMESCF*38, FNAMEDEL*38, FNAMEUP*38, FNAMEDOWN*38,
     .  PASTE*38

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE, PLANE,
     .  NEIGHB, WROUT

C **********************************************************************
C INTEGER NPLAMAX          : Maximum number of points in the plane
C REAL*8  PLAPO(NPLAMAX,3) : Coordinates of the points of the plane in PRF
C REAL*8  POINRE(NPLAMAX,3): Coordinates of the points of the plane in Lattice
C                            Reference Frame
C INTEGER NAPLA            : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA)      : Indices of the atoms whose coordinates will 
C                            be rotated from the lattice reference frame 
C                            to the in-plane reference frame
C REAL*8  XAPLA(3,NA)      : Atomic coordinates in plane reference frame
C INTEGER IA               : Atom whose neighbours are needed.
C                            A routine initialization must be done by
C                            a first call with IA = 0
C                            If IA0=0, point X0 is used as origin instead
C INTEGER ISEL             : Single-counting switch (0=No, 1=Yes). If ISEL=1,
C                            only neighbours with JA.LE.IA are included in JNA
C INTEGER NNA              : Number of non-zero orbitals at a point in 
C                            real space
C INTEGER JNA(MAXNA)       : Atom index of neighbours. The neighbours
C                            atoms might be in the supercell
C REAL*8  XIJ(3,MAXNA)     : Vectors from point in real space to orbitals
C REAL*8  R2IJ(MAXNA)      : Squared distance to atomic orbitals
C REAL*8  XPO(3)           : Coordinates of the point of the plane respect
C                            we are going to calculate the neighbours orbitals
C REAL*8  DENCHAR          : Self-Consistent Density of Charge at a given
C                            point in real space
C REAL*8  DENCHAR0         : Harris Density of Charge at a given point
C                            in real space
C **********************************************************************

C Allocate some variables ---------------------------------------------
      NPLAMAX = NPX * NPY

      ALLOCATE(PLAPO(NPLAMAX,3))
      CALL MEMORY('A','D',3*NPLAMAX,'rhoofr')

      ALLOCATE(POINRE(NPLAMAX,3))
      CALL MEMORY('A','D',3*NPLAMAX,'rhoofr')

      ALLOCATE(INDICES(NA))
      CALL MEMORY('A','I',NA,'rhoofr')

      ALLOCATE(XAPLA(3,NA))
      CALL MEMORY('A','D',3*NA,'rhoofr')

      ALLOCATE(DISCF(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

      ALLOCATE(DIATM(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

      ALLOCATE(DIUP(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

      ALLOCATE(DIDOWN(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

C Build the plane ------------------------------------------------------
          CALL PLANE( NA, NPLAMAX, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                XA, NAPLA, INDICES, ISCALE,
     .                POINRE, PLAPO, XAPLA )   
C End building of the plane --------------------------------------------

C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = MAXNA
      IF (ALLOCATED(JNA)) THEN
        CALL MEMORY('D','I',SIZE(JNA),'rhoofr')
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        CALL MEMORY('D','D',SIZE(R2IJ),'rhoofr')
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        CALL MEMORY('D','D',SIZE(XIJ),'rhoofr')
        DEALLOCATE(XIJ)
      ENDIF
      ALLOCATE(JNA(MAXNA))
      CALL MEMORY('A','I',MAXNA,'rhoofr')
      ALLOCATE(R2IJ(MAXNA))
      CALL MEMORY('A','D',MAXNA,'rhoofr')
      ALLOCATE(XIJ(3,MAXNA))
      CALL MEMORY('A','D',3*MAXNA,'rhoofr')

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

C Open files to store charge density -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')
      IF (NSPIN .EQ. 1) THEN
        FNAMESCF = PASTE(SNAME,'.CON.SCF')
        FNAMEDEL = PASTE(SNAME,'.CON.DEL')
        CALL IO_ASSIGN(UNIT1)
        OPEN(UNIT = UNIT1, FILE = FNAMESCF, STATUS = 'UNKNOWN',
     .       FORM = 'FORMATTED')
        REWIND(UNIT1)
        CALL IO_ASSIGN(UNIT2)
        OPEN(UNIT = UNIT2, FILE = FNAMEDEL, STATUS = 'UNKNOWN',
     .       FORM = 'FORMATTED')
        REWIND(UNIT2)
      ELSEIF (NSPIN .EQ. 2) THEN
        FNAMESCF = PASTE(SNAME,'.CON.MAG' )
        FNAMEDEL = PASTE(SNAME,'.CON.DEL' )
        FNAMEUP  = PASTE(SNAME,'.CON.UP'  )
        FNAMEDOWN= PASTE(SNAME,'.CON.DOWN')
        CALL IO_ASSIGN(UNIT1)
        OPEN(UNIT = UNIT1, FILE = FNAMESCF, STATUS = 'UNKNOWN',
     .       FORM = 'FORMATTED')
        REWIND(UNIT1)
        CALL IO_ASSIGN(UNIT2)
        OPEN(UNIT = UNIT2, FILE = FNAMEDEL, STATUS = 'UNKNOWN',
     .       FORM = 'FORMATTED')
        REWIND(UNIT2)
        CALL IO_ASSIGN(UNIT3)
        OPEN(UNIT = UNIT3, FILE = FNAMEUP, STATUS = 'UNKNOWN',
     .       FORM = 'FORMATTED')
        REWIND(UNIT3)
        CALL IO_ASSIGN(UNIT4)
        OPEN(UNIT = UNIT4, FILE = FNAMEDOWN, STATUS = 'UNKNOWN',
     .       FORM = 'FORMATTED')
        REWIND(UNIT4)
      ELSE
        WRITE(6,*)'BAD NUMBER NSPIN IN RHOOFR.F'
        WRITE(6,*)'NSPIN = ',NSPIN
        WRITE(6,*)'IT MUST BE 1 => NON POLARIZED, OR 2 = > POLARIZED'
        STOP
      ENDIF

   
C Dump input data in output files --------------------------------------
      CALL WROUT( IOPTION, UNIT1, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .            NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD,
     .            NA, NAPLA, INDICES, XAPLA )
      CALL WROUT( IOPTION, UNIT2, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .            NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD, 
     .            NA, NAPLA, INDICES, XAPLA )
      IF ( NSPIN .EQ. 2) THEN
        CALL WROUT( IOPTION, UNIT3, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .              NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD, 
     .              NA, NAPLA, INDICES, XAPLA )
        CALL WROUT( IOPTION, UNIT4, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .              NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD, 
     .              NA, NAPLA, INDICES, XAPLA )
      ENDIF
      
C Loop over all points in real space -----------------------------------
      DO 100 NPO = 1, NPX*NPY
C Initialize the density of charge at each point -----------------------
        DENCHAR  = 0.D0
        DENCHAR0 = 0.D0
        DENUP    = 0.D0
        DENDOWN  = 0.D0

C Localize non-zero orbitals at each point in real space ---------------
        DO IX = 1,3
          XPO(IX) = POINRE(NPO,IX)
        ENDDO
     
        IA   = 0
        ISEL = 0
        NNA  = MAXNA

        CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .               NNA, JNA, XIJ, R2IJ, FIRST )

C Loop over Non-zero orbitals ------------------------------------------ 
         DO 110 IAT1 = 1, NNA
           IF( R2IJ(IAT1) .GT. RMAX2 ) GOTO 110

           IAVEC1   = JNA(IAT1)
           IS1      = ISA(IAVEC1)
           XVEC1(1) = -XIJ(1,IAT1)
           XVEC1(2) = -XIJ(2,IAT1)
           XVEC1(3) = -XIJ(3,IAT1)


           DO 120 IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
             IPHI1 = IPHORB(IO)
             IUO   = INDXUO(IO)
             CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

C Copy full row IUO of Density Matrix to DI(j) --------------------------
             IF (IO . EQ. IUO) THEN
               DO 130 IN = 1, NUMD(IUO)
                 IND = IN + LISTDPTR(IUO)
                 J = LISTD(IND)
                 IF (NSPIN .EQ. 1) THEN
                   DISCF(J) = DSCF(IND,1)
                   DIATM(J) = DSCFNA(IND)
                 ELSEIF (NSPIN .EQ. 2) THEN
                   DIUP(J)   = DSCF(IND,1)
                   DIDOWN(J) = DSCF(IND,2)
                   DIATM(J)  = DSCFNA(IND)
                 ENDIF
 130           ENDDO
             ELSE
               DO 135 IN = 1, NUMD(IUO)
                 IND = IN + LISTDPTR(IUO)
                 J = LISTSC( IO, IUO, LISTD(IND) )
                 IF (NSPIN .EQ. 1) THEN
                   DISCF(J) = DSCF(IND,1)
                   DIATM(J) = DSCFNA(IND)
                 ELSEIF (NSPIN .EQ. 2) THEN
                   DIUP(J)   = DSCF(IND,1)
                   DIDOWN(J) = DSCF(IND,2)
                   DIATM(J)  = DSCFNA(IND)
                 ENDIF
 135           ENDDO
             ENDIF

C Loop again over neighbours -------------------------------------------
             DO 140 IAT2 = 1, NNA
               IAVEC2   = JNA(IAT2)
               IS2      = ISA(IAVEC2)
               XVEC2(1) = -XIJ(1,IAT2)
               XVEC2(2) = -XIJ(2,IAT2)
               XVEC2(3) = -XIJ(3,IAT2)
               
               DO 150 JO = LASTO(IAVEC2-1) + 1, LASTO(IAVEC2)
                 IPHI2 = IPHORB(JO)
                 CALL PHIATM( IS2, IPHI2, XVEC2, PHINU, GRPHINU )
                 
                 IF ( NSPIN .EQ. 1 ) THEN
                   DENCHAR  = DENCHAR  + PHINU*PHIMU*DISCF(JO)
                   DENCHAR0 = DENCHAR0 + PHINU*PHIMU*DIATM(JO)
                 ELSEIF (NSPIN .EQ. 2) THEN 
                   DENUP    = DENUP    + PHINU*PHIMU*DIUP(JO)
                   DENDOWN  = DENDOWN  + PHINU*PHIMU*DIDOWN(JO)
                   DENCHAR0 = DENCHAR0 + PHINU*PHIMU*DIATM(JO)
                 ENDIF
                 
 150           ENDDO

 140         ENDDO
 120       ENDDO
 110     ENDDO
 
         IF ( NSPIN .EQ. 1 ) THEN
           WRITE(UNIT1,'(3F12.5)')
     .          PLAPO(NPO,1),PLAPO(NPO,2),DENCHAR*ARMUNI
           WRITE(UNIT2,'(3F12.5)')
     .          PLAPO(NPO,1),PLAPO(NPO,2),(DENCHAR-DENCHAR0)*ARMUNI 
         ELSEIF ( NSPIN .EQ. 2 ) THEN
           WRITE(UNIT1,'(3F12.5)')
     .          PLAPO(NPO,1),PLAPO(NPO,2),(DENUP-DENDOWN)*ARMUNI
           WRITE(UNIT2,'(3F12.5)')
     .          PLAPO(NPO,1),PLAPO(NPO,2),
     .          (DENUP+DENDOWN-DENCHAR0)*ARMUNI
           WRITE(UNIT3,'(3F12.5)')
     .          PLAPO(NPO,1),PLAPO(NPO,2),DENUP*ARMUNI
           WRITE(UNIT4,'(3F12.5)')
     .          PLAPO(NPO,1),PLAPO(NPO,2),DENDOWN*ARMUNI
         ENDIF

         IF ( MOD(NPO,NPX) .EQ. 0 ) THEN
           WRITE(UNIT1,'(\)')
           WRITE(UNIT2,'(\)')
           IF ( NSPIN .EQ. 2 ) THEN
             WRITE(UNIT3,'(\)')
             WRITE(UNIT4,'(\)')
           ENDIF
         ENDIF

 100  ENDDO  

      IF ( NSPIN .EQ. 1) THEN
        CALL IO_CLOSE(UNIT1)
        CALL IO_CLOSE(UNIT2)
      ELSEIF (NSPIN .EQ. 2) THEN
        CALL IO_CLOSE(UNIT1)
        CALL IO_CLOSE(UNIT2)
        CALL IO_CLOSE(UNIT3)
        CALL IO_CLOSE(UNIT4)
      ENDIF
     
          
      RETURN    
      END
