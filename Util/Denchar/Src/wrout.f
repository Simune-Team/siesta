      SUBROUTINE WROUT(IOPTION, UNIT1, NORMAL, COORPO, DIRVER1, DIRVER2, 
     .                 NPX, NPY, XMIN, XMAX, YMIN, YMAX, IUNITCD,
     .                 MAXATOM, NAPLA, INDICES, XAPLA )
C **********************************************************************
C Dump input data into the output files
C Written by J. Junquera Feb '99
C **********************************************************************

      IMPLICIT NONE

      INCLUDE 'fdfdefs.h'

      INTEGER
     .  IOPTION, UNIT1, NPX, NPY, IUNITCD, MAXATOM, NAPLA, 
     .  INDICES(MAXATOM)

      DOUBLE PRECISION
     .  NORMAL(3), COORPO(3,3), DIRVER1(3), DIRVER2(3), XAPLA(3,MAXATOM)

      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX

C **************  INPUT  ***********************************************
C INTEGER IOPTION        : Option to generate the plane
C                          1 = Normal vector
C                          2 = Two vectors belonging to the plane
C                          3 = Three points of the plane
C                          4 = Three atomic indices
C INTEGER UNIT1          : Number of the logical unit where the data will 
C                          be dumped
C REAL*8  NORMAL(3)      : Components of the normal vector
C REAL*8  COORPO(3,3)    : Coordinates of the three points used to define
C                          the plane 
C REAL*8  DIRVER(3)      : Components of two vector contained in the plane
C INTEGER NPX,NPY        : Number of points generated along x and y
C                          direction in a system of reference in which
C                          the third components od the points of the plane is
C                          zero (Plane Reference Frame; PRF)
C REAL*8  XMIN, XMAX     : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX     : Limits of the plane in the PRF for y-direction
C INTEGER IUNITCD        : Unit of the charge density
C INTEGER MAXATOM        : Maximum number of atoms
C INTEGER NAPLA          : Number of atoms whose coordiantes has been rotated   
C INTEGER INDICES(MAXATOM): Indices of tha atoms whose coordinates has 
C                           been roated
C REAL*8  XAPLA(3,MAXATOM): Atomic coordiantes in the in-plane reference frame
C **********************************************************************

C ***************  INTERNAL VARIABLES **********************************
      CHARACTER*30
     .  SNAME

      INTEGER
     .  IX, IP, IA

C Open files to store charge density -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')

      WRITE(UNIT1,'(A)')
     .    '#                          ************************       '
      WRITE(UNIT1,'(A)')
     .    '#                          *  WELCOME TO DENCHAR  *       '
      WRITE(UNIT1,'(A)')
     .    '#                          ************************       '

      WRITE(UNIT1,'(A,A)')
     .    '#  WROUT: You are running DENCHAR for system: ',SNAME
      WRITE(UNIT1,'(A,/A)')
     .    '#  WROUT: The options you have chosen to generate the plane',
     .    '#  WROUT: are the following: '

      IF( IOPTION .EQ. 1 ) THEN

        WRITE(UNIT1,'(A)')
     .    '#  WROUT: Option to generate the plane : NormalVector'
        WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '#  WROUT: Components of the normal vector : ',
     .    '#  WROUT: ',(NORMAL(IX),IX=1,3)
        WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '#  WROUT: Origin of the plane : ',
     .    '#  WROUT: ',(COORPO(1,IX),IX=1,3)
        WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '#  WROUT: Another point to define the X direction : ',
     .    '#  WROUT: ',(COORPO(2,IX),IX=1,3)

      ELSEIF( IOPTION .EQ. 2 ) THEN 
        
        WRITE(UNIT1,'(A)')
     .    '#  WROUT: Option to generate the plane : TwoLines'
        WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '#  WROUT: Components of the first vector inside the plane :',
     .    '#  WROUT: ',(DIRVER1(IX),IX=1,3)
        WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '#  WROUT: Components of the second vector inside the plane:',
     .    '#  WROUT: ',(DIRVER2(IX),IX=1,3)
        WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '#  WROUT: Origin of the plane : ',
     .    '#  WROUT: ',(COORPO(1,IX),IX=1,3)

      ELSEIF( IOPTION .EQ. 3 ) THEN 

        WRITE(UNIT1,'(A)')
     .    '#  WROUT: Option to generate the plane : ThreePoints'
        WRITE(UNIT1,'(A)')
     .    '#  WROUT: Coordinates of three points in the plane : '
        DO IP = 1,3
          WRITE(UNIT1,'(A,3F12.5)')
     .      '#  WROUT:',(COORPO(IP,IX),IX=1,3)
        ENDDO
      
      ELSEIF( IOPTION .EQ. 4 ) THEN 

        WRITE(UNIT1,'(A)')
     .    '#  WROUT: Option to generate the plane : ThreeAtomicIndices'
        WRITE(UNIT1,'(A)')
     .    '#  WROUT: Position of the three atoms : '
        DO IP = 1,3
          WRITE(UNIT1,'(A,3F12.5)')
     .      '#  WROUT:',(COORPO(IP,IX),IX=1,3)
        ENDDO

      ENDIF

      WRITE(UNIT1,'(A,/,A,I5)')
     .    '#  WROUT: Number of points in the x-direction : ',
     .    '#  WROUT: ', NPX
      WRITE(UNIT1,'(A,/,A,I5)')
     .    '#  WROUT: Number of points in the y-direction : ',
     .    '#  WROUT: ', NPY
      WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '#  WROUT: Minimum value of the x-component of the window : ',
     .    '#  WROUT: ', XMIN,' bohrs'
      WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '#  WROUT: Maximum value of the x-component of the window : ',
     .    '#  WROUT: ', XMAX,' bohrs'
      WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '#  WROUT: Minimum value of the y-component of the window : ',
     .    '#  WROUT: ', YMIN,' bohrs'
      WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '#  WROUT: Maximum value of the y-component of the window : ',
     .    '#  WROUT: ', YMAX,' bohrs'

      IF ( IUNITCD .EQ. 1) THEN
        WRITE(UNIT1,'(A,/,A)')
     .    '#  WROUT: Unit of the charge density in output files : ',
     .    '#  WROUT: Electrons/(bohr**3)'
      ELSEIF ( IUNITCD .EQ. 2) THEN
        WRITE(UNIT1,'(A,/,A)')
     .    '#  WROUT: Unit of the charge density in output files : ',
     .    '#  WROUT: Electrons/(angstrom**3)'
      ELSEIF( IUNITCD .EQ. 3) THEN
        WRITE(UNIT1,'(A,/,A)')
     .    '#  WROUT: Unit of the charge density in output files : ',
     .    '#  WROUT: Electrons/unit cell'
      ENDIF

      IF( NAPLA .NE. 0) THEN
        WRITE(UNIT1,'(A)')
     .    '#  WROUT: Atomic coordinates in the in-plane reference frame'
        WRITE(UNIT1,'(A,19(1H ),A)')
     .    '#  WROUT: Atomic Index','Atomic coordinates'
        DO IA = 1, NAPLA
          WRITE(UNIT1,'(A,I14,5X,3F15.4)')
     .      '#',INDICES(IA), (XAPLA(IX,INDICES(IA)),IX=1,3)
        ENDDO
      ENDIF

      WRITE(UNIT1,'(A)')

      END
       
      
