
      SUBROUTINE RHOOFR( CELL, NA, NAMAX, MAXNO, MAXA, MAXO, MAXSPN,
     .                   NSPIN, RMAXO, LASTO, XA, ISA, IPHORB,
     .                   INDXUO, NUMD, LISTD, DSCF, DSCFNA, 
     .                   IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                   NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                   ARMUNI, IUNITCD )
C **********************************************************************
C Compute the density of charge at the points of a plane in real space
C Coded by J. Junquera November'98
C **********************************************************************

      IMPLICIT NONE

      INCLUDE 'fdfdefs.h'

C INTEGER ORBMAX    : Maximum number of Atomic Orbitals
C INTEGER ATOMAX    : Maximum number of neighbour atoms
C INTEGER NPLAMAX   : Maximum number of points in the plane

      INTEGER NPLAMAX, ORBMAX, ATOMAX
      PARAMETER(ORBMAX  = 10000 )
      PARAMETER(ATOMAX  =  1000 )
      PARAMETER(NPLAMAX = 10000 )

      INTEGER
     .  IOPTION, NPX, NPY, MAXA, MAXO,
     .  NAMAX, MAXNO, MAXSPN, NA, NSPIN, IUNITCD,
     .  LASTO(0:MAXA), ISA(MAXA), IPHORB(MAXO),
     .  INDXUO(MAXO), NUMD(MAXO), LISTD(MAXNO,MAXO)

      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), ARMUNI,
     .  PLAPO(NPLAMAX,3), POINRE(NPLAMAX,3)

      DOUBLE PRECISION
     . CELL(3,3), XA(3,MAXA), RMAXO, DSCF(MAXNO,MAXO,MAXSPN),
     . DSCFNA(MAXNO,MAXO,MAXSPN)

      DOUBLE PRECISION
     .  DENCHAR, DENCHAR0, DENUP, DENDOWN
 
      LOGICAL FIRST

C **********************************************************************
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
C REAL*8  DIRVER(3)        : Components of two vector contained in the plane
C                            (Only if ioption = 2)
C REAL*8  COORPO(3,3)      : Coordinates of the three points used to define
C                            the plane (Only if ioption = 3)
C REAL*8  PLAPO(NPLAMAX,3) : Coordinates of the points of the plane in PRF
C REAL*8  POINRE(NPLAMAX,3): Coordinates of the points of the plane in Lattice
C                            Reference Frame
C INTEGER NAMAX            : Maximum number of neighbour atoms of any atom
C INTEGER MAXNO            : Maximum number of neighbours orbitals 
C INTEGER MAXA             : Maximum number of atoms
C INTEGER MAXO             : Maximum number of atomic orbitals
C INTEGER MAXSPN           : Maximum number of spin
C INTEGER ISA(MAXA)        : Species index of each atom
C INTEGER IPHORB(MAXO)     : Orital index of each orbital in its atom
C INTEGER INDXUO(MAXO)     : Equivalent orbital in unit cell
C INTEGER NA               : Total number of atoms in Supercell
C INTEGER NSPIN            : Number of different spin polarizations
C                            Nspin = 1 => unpolarized, Nspin = 2 => polarized
C INTEGER LASTO(0:MAXA)    : Last orbital of each atom in array iphorb
C INTEGER NUMD(MAXO)       : Control vector of Density Matrix
C                            (number of non-zero elements of eaach row)
C INTEGER LISTD(MAXNO,MAXO): Control vector of Density Matrix
C                            (list of non-zero elements of each row)
C INTEGER IUNITCD          : Unit of the charge density
C REAL*8  CELL(3,3)        : Unit cell vectors CELL(IXYZ,IVECT)
C REAL*8  XA(3,NA)         : Atomic positions in cartesian coordinates
C REAL*8  RMAXO            : Maximum range of basis orbitals
C REAL*8  DSCF(MAXNO,MAXO,MAXPSN)  : Density Matrix
C REAL*8  DSCFNA(MAXNO,MAXO,MAXPSN): Density Matrix for Neutral Atoms
C REAL*8  DENCHAR          : Self-Consistent Density of Charge at a given
C                            point in real space
C REAL*8  DENCHAR0         : Harris Density of Charge at a given point
C                            in real space
C **********************************************************************

      INTEGER
     .  NPO, IA, ISEL, NNA, JANA(ATOMAX), UNIT1, UNIT2, UNIT3, UNIT4
 
      INTEGER
     .  I, J, IX, IN, IAT1, IAT2, JO, IO, IUO, IAVEC, IAVEC1, IAVEC2, 
     .  IS1, IS2, IPHI1, IPHI2

      DOUBLE PRECISION
     .  RMAX, XIJ(3,ATOMAX), R2IJ(ATOMAX), XPO(3), RMAX2, XVEC1(3),
     .  XVEC2(3), PHIMU, PHINU, GRPHIMU(3), GRPHINU(3), DISCF(ORBMAX),
     .  DIATM(ORBMAX), DIUP(ORBMAX), DIDOWN(ORBMAX)

      CHARACTER
     .  SNAME*30, FNAMESCF*38, FNAMEDEL*38, FNAMEUP*38, FNAMEDOWN*38,
     .  PASTE*38
 

      EXTERNAL
     .  CHKDIM, IO_ASSIGN, IO_CLOSE, PASTE, PLANE,
     .  NEIGHB, PHIATM, WROUT

C **********************************************************************
C INTEGER IA               : Atom whose neighbours are needed.
C                            A routine initialization must be done by
C                            a first call with IA = 0
C INTEGER ISEL             : Single-counting switch (0=No, 1=Yes). If ISEL=1,
C                            only neighbours with JA.LE.IA are included in JANA
C INTEGER NNA              : Number of non-zero orbitals at a point in 
C                            real space
C INTEGER JANA(NAMAX)      : Atom index of neighbours
C REAL*8  XIJ(3,NAMAX)     : Vectors from point in real space to orbitals
C REAL*8  R2IJ(NAMAX)      : Squared distance to atomic orbitals
C REAL*8  XPO(3)           : Coordinates of the point of the plane respect
C                            we are going to calculate the neighbours orbitals
C **********************************************************************

C Check some dimensions ------------------------------------------------
          CALL CHKDIM( 'RHOOFR', 'NPLAMAX', NPLAMAX, NPX*NPY, 1 )
          CALL CHKDIM( 'RHOOFR', 'ATOMAX' , ATOMAX , NAMAX  , 1 )
          CALL CHKDIM( 'RHOOFR', 'ORBMAX' , ORBMAX , MAXO   , 1 )

C Build the plane ------------------------------------------------------
          CALL PLANE( NPLAMAX, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                POINRE, PLAPO )   
C End building of the plane --------------------------------------------

C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = NAMAX
      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JANA, XIJ, R2IJ, FIRST )
 
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
     .            NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD )
      CALL WROUT( IOPTION, UNIT2, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .            NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD )
      IF ( NSPIN .EQ. 2) THEN
        CALL WROUT( IOPTION, UNIT3, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .              NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD )
        CALL WROUT( IOPTION, UNIT4, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .              NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD )
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
        NNA  = NAMAX

        CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .               NNA, JANA, XIJ, R2IJ, FIRST )

C Loop over Non-zero orbitals ------------------------------------------ 
         DO 110 IAT1 = 1, NNA
           IF( R2IJ(IAT1) .GT. RMAX2 ) GOTO 110

           IAVEC1   = JANA(IAT1)
           IS1      = ISA(IAVEC1)
           XVEC1(1) = -XIJ(1,IAT1)
           XVEC1(2) = -XIJ(2,IAT1)
           XVEC1(3) = -XIJ(3,IAT1)

           DO 120 IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
             IPHI1 = IPHORB(IO)
             CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

C Copy full row IUO of Density Matrix to DI(j) --------------------------
             DO 130 IN = 1, NUMD(IO)
               J       = LISTD(IN,IO)
               IUO     = INDXUO(IO)
               IF (NSPIN .EQ. 1) THEN
                 DISCF(J) = DSCF(IN,IUO,1)
                 DIATM(J) = DSCFNA(IN,IUO,1) 
               ELSEIF (NSPIN .EQ. 2) THEN
                 DIUP(J)   = DSCF(IN,IUO,1)
                 DIDOWN(J) = DSCF(IN,IUO,2)
                 DIATM(J)  = DSCFNA(IN,IUO,1)
               ENDIF
 130         ENDDO
C Loop again over neighbours -------------------------------------------
             DO 140 IAT2 = 1, NNA
               IAVEC2   = JANA(IAT2)
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
