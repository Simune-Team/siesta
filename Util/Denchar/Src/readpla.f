
      SUBROUTINE READPLA( IOPTION, XMIN, XMAX, YMIN, YMAX,
     .                    NPX, NPY, COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                    ARMUNI, MAXA, XA, VOLUME, IUNITCD, ISCALE ) 
C **********************************************************************
C Read the data file to prepare the plane in which we are going to
C calculate the charge density
C Coded by J. Junquera, November 98
C **********************************************************************

      IMPLICIT NONE

      INCLUDE 'fdfdefs.h'

      INTEGER 
     .  IOPTION, NPX, NPY, ISCALE, IUNITCD
     
      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX, 
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3),
     .  ARMUNI

      INTEGER
     .  MAXA
      
      DOUBLE PRECISION
     .  XA(3,MAXA), VOLUME

C **********************************************************************
C INTEGER IOPTION        : Option to generate the plane
C                          1 = Normal vector
C                          2 = Two vectors belonging to the plane
C                          3 = Three points of the plane
C                          4 = Three atomic indices
C INTEGER NPX, NPY       : Number of points generated along x and y
C                          directions ina a system of reference in which
C                          the third component of the points of the plane
C                          is zero (Plane Reference Frame; PRF)
C REAL*8  XMIN, XMAX     : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX     : Limits of the plane in the PRF for y-direction
C REAL*8  NORMAL(3)      : Components of the normal vector used to define
C                          the plane 
C REAL*8  DIRVER(3)      : Components of the two vectors contained in the plane 
C REAL*8  COORPO(3,3)    : Coordinates of the three points used to define 
C                          the plane. COORPO(POINT,IX)
C INTEGER MAXA           : Maximum number of atoms
C REAL*8  XA(3,MAXA)     : Atomic coordinates
C REAL*8 VOLUME          : Volumen of unit cell (in bohr**3)
C INTEGER IUNITCD        : Units for the electron density
C                          IUNITCD = 1 => Ele/(bohr)**3
C                          IUNITCD = 2 => Ele/(Ang)**3
C                          IUNITCD = 3 => Ele/(unitcell)
C INTEGER ISCALE         : Units for the atomic positions
C                          (ISCALE = 1 => Bohrs, ISCALE = 2 => Ang)
C **********************************************************************

      CHARACTER 
     .  OGP*22, OGP_DEFECT*22,
     .  CPF*22, CPF_DEFECT*22,
     .  UCD*22, UCD_DEFECT*22

      INTEGER
     .  IUNIT, IX, JX, NPX_DEFECT, NPY_DEFECT, IND1, IND2, IND3

      DOUBLE PRECISION
     .  ORIGIN(3), XDIR(3)

      LOGICAL 
     .  LEQI, COLIN

      EXTERNAL 
     .  LEQI, COLINEAR

      DATA ORIGIN /0.D0,0.D0,0.D0/
      DATA XDIR   /1.D0,0.D0,0.D0/
      DATA IND1   /1/
      DATA IND2   /2/
      DATA IND3   /3/
      DATA COLIN  /.FALSE./

      CPF_DEFECT = 'Bohr'
      CPF = FDF_STRING('2D.CoorUnits',CPF_DEFECT)
      IF (LEQI(CPF,'bohr')) then
        ISCALE = 1
      ELSEIF (LEQI(CPF,'ang')) then
        ISCALE = 2
      ENDIF

      UCD_DEFECT = 'Ele/bohr**3'
      UCD = FDF_STRING('2D.DensityUnits',UCD_DEFECT)
      IF (LEQI(UCD,'ele/bohr**3')) then
        IUNITCD = 1
      ELSEIF (LEQI(UCD,'ele/ang**3')) then
        IUNITCD = 2
      ELSEIF (LEQI(UCD,'ele/unitcell')) then
        IUNITCD = 3
      ELSE
       WRITE(6,'(A)')' readpla:  Wrong Option in Units of             '
       WRITE(6,'(A)')' readpla:        Charge Density                 '
       WRITE(6,'(A)')' readpla:  You must choose one of the following:'       
       WRITE(6,'(A)')'                                                '
       WRITE(6,'(A)')' readpla:      - Ele/bohr**3                    '
       WRITE(6,'(A)')' readpla:      - Ele/ang**3                     '
       WRITE(6,'(A)')' readpla:      - Ele/unitcell                   '
       STOP
      ENDIF

      NPX_DEFECT = 50
      NPY_DEFECT = 50
      NPX = FDF_INTEGER('2D.NumberPointsX',NPX_DEFECT)
      NPY = FDF_INTEGER('2D.NumberPointsY',NPY_DEFECT)

      XMIN = FDF_PHYSICAL('2D.MinX',-3.D0,'Bohr')
      XMAX = FDF_PHYSICAL('2D.MaxX', 3.D0,'Bohr')
      YMIN = FDF_PHYSICAL('2D.MinY',-3.D0,'Bohr')
      YMAX = FDF_PHYSICAL('2D.MaxY', 3.D0,'Bohr')

      OGP_DEFECT = 'NormalVector'
      OGP = FDF_STRING('2D.PlaneGeneration',OGP_DEFECT)
      IF (LEQI(OGP,'normalvector')) then
        IOPTION = 1
      ELSEIF (LEQI(OGP,'twolines')) then
        IOPTION = 2
      ELSEIF (LEQI(OGP,'threepoints')) then
        IOPTION = 3
      ELSEIF (LEQI(OGP,'threeatomicindices')) then
        IOPTION = 4
      ELSE
       WRITE(6,'(A)')' readpla:  Wrong Option to Generate The Plane   '
       WRITE(6,'(A)')' readpla:  You must choose one of the following:'       
       WRITE(6,'(A)')'                                                '
       WRITE(6,'(A)')' readpla:      - NormalVector                   '
       WRITE(6,'(A)')' readpla:      - TwoLines                       '
       WRITE(6,'(A)')' readpla:      - ThreePoints                    '
       WRITE(6,'(A)')' readpla:      - ThreeAtomicIndices             '
       STOP
      ENDIF

      IF ( FDF_BLOCK('2D.CompNormalVector',IUNIT) ) THEN
        READ(IUNIT,*)(NORMAL(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('2D.Comp2Vectors',IUNIT) ) THEN
        READ(IUNIT,*)(DIRVER1(IX),IX=1,3)
        READ(IUNIT,*)(DIRVER2(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('2D.Coor3Points',IUNIT) ) THEN
        READ(IUNIT,*)(COORPO(1,IX),IX=1,3)
        READ(IUNIT,*)(COORPO(2,IX),IX=1,3)
        READ(IUNIT,*)(COORPO(3,IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('2D.Indices3Atoms',IUNIT) ) THEN
        READ(IUNIT,*)IND1, IND2, IND3
      ENDIF

      IF ( IOPTION .EQ. 4 ) THEN
        DO IX = 1,3
          COORPO(1,IX) = XA(IX,IND1)
        ENDDO
        DO IX = 1,3
          COORPO(2,IX) = XA(IX,IND2)
        ENDDO
        DO IX = 1,3
          COORPO(3,IX) = XA(IX,IND3)
        ENDDO
      ENDIF

C Check if the three points are colinear -------------------------------
      IF ((IOPTION .EQ. 3) .OR. (IOPTION .EQ. 4))THEN
         CALL COLINEAR( COORPO, COLIN )
         IF(COLIN) THEN
           WRITE(6,*)'The coordinates of the three points are colinear'
           WRITE(6,*)'and do not define a plane' 
           WRITE(6,*)'Please, check these coordinates in the input file'      
           STOP
         ENDIF
      ENDIF

      IF ( FDF_BLOCK('2D.PlaneOrigin',IUNIT) ) THEN
        READ(IUNIT,*)(ORIGIN(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('2D.X_Axis',IUNIT) ) THEN
        READ(IUNIT,*)(XDIR(IX),IX=1,3)
      ENDIF

      IF (IOPTION .LT. 3) THEN
        DO IX = 1,3      
          COORPO(1,IX) = ORIGIN(IX)
        ENDDO
        IF(IOPTION .EQ. 1) THEN
          DO IX = 1,3
            COORPO(2,IX) = XDIR(IX)
          ENDDO
        ENDIF
      ENDIF

C Scale points coordinates
C   Iscale = 1 => Do nothing
C   Iscale = 2 => Multiply by 1./0.529177 (Ang --> Bohr)

      IF( (ISCALE .EQ. 2) .AND. (IOPTION .NE. 4) ) THEN
        DO IX = 1,3
          DO JX = 1,3
            COORPO(JX,IX) = 1.D0 / 0.529177D0 * COORPO(JX,IX)
          ENDDO
          ORIGIN(IX)  = 1.D0 / 0.529177D0 * ORIGIN(IX)
          DIRVER1(IX) = 1.D0 / 0.529177D0 * DIRVER1(IX)
          DIRVER2(IX) = 1.D0 / 0.529177D0 * DIRVER2(IX)
        ENDDO
      ENDIF 

C Units of Charge Density
C   Iunitcd = 1 => Do nothing
C   Iunitcd = 2 => Multiply by (1.d0 / 0.529177d0) **3 (bohr**3 --> Ang**3)
C   Iunitcd = 3 => Multiply by volume unit cell (in bohrs**3) 

      IF (IUNITCD .EQ. 1) THEN
        ARMUNI = 1.D0
      ELSEIF( IUNITCD .EQ. 2 ) THEN
        ARMUNI = (1.D0 / 0.529177D0)**3 
      ELSEIF( IUNITCD .EQ. 3 ) THEN
        ARMUNI = VOLUME
      ENDIF

      END

