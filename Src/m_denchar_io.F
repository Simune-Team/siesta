      module m_denchar_io

      use m_denchar_geom
      use ionew

      public :: readpla

      CONTAINS


      SUBROUTINE READPLA(MAXA, XA, VOLUME, IDIMEN,
     .                  IOPTION, IUNITCD, ISCALE, NPX, NPY, NPZ,
     .                  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .                  COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                  ARMUNI ) 

C **********************************************************************
C Read the data file to prepare the plane or 3D grid 
C in which we are going to  calculate the charge density
C or wavefunctions
C
C Coded by J. Junquera, November 98
C Modified by P. Ordejon to include 3D and wavefunction capabilities
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  MAXA, IDIMEN
      
      DOUBLE PRECISION, INTENT(IN) ::
     .  XA(3,MAXA), VOLUME

      INTEGER, INTENT(OUT) ::
     .  IOPTION, NPX, NPY, NPZ, ISCALE, IUNITCD
     
      DOUBLE PRECISION, INTENT(OUT) ::
     .  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), 
     .  ARMUNI

C ****** INPUT *********************************************************
C INTEGER MAXA           : Maximum number of atoms
C REAL*8  XA(3,MAXA)     : Atomic coordinates
C REAL*8 VOLUME          : Volumen of unit cell (in bohr**3)
C INTEGER IDIMEN         : Specify if the run is to plot quantities
C                          in a plane or in a 3D grid (2 or 3, respect)
C ****** OUTPUT ********************************************************
C INTEGER IOPTION        : Option to generate the plane
C                          1 = Normal vector
C                          2 = Two vectors belonging to the plane
C                          3 = Three points of the plane
C                          4 = Three atomic indices
C INTEGER IUNITCD        : Units for the electron density
C                          IUNITCD = 1 => Ele/(bohr)**3
C                          IUNITCD = 2 => Ele/(Ang)**3
C                          IUNITCD = 3 => Ele/(unitcell)
C INTEGER ISCALE         : Units for the atomic positions
C                          (ISCALE = 1 => Bohrs, ISCALE = 2 => Ang)
C INTEGER NPX, NPY, NPZ  : Number of points generated along x and y
C                          (and z, for 3D-grids) directions ina a system 
C                          of reference in which the third component of 
C                          the points of the plane is zero (Plane 
C                          Reference Frame; PRF)
C REAL*8  XMIN, XMAX     : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX     : Limits of the plane in the PRF for y-direction
C REAL*8  ZMIN, ZMAX     : Limits of the z-direction in the PRF (for 3D-grids)
C REAL*8  COORPO(3,3)    : Coordinates of the three points used to define 
C                          the plane. COORPO(POINT,IX)
C REAL*8  NORMAL(3)      : Components of the normal vector used to define
C                          the plane 
C REAL*8  DIRVER1(3)     : Components of the first vector 
C                          contained in the xy plane 
C REAL*8  DIRVER2(3)     : Components of the first vector 
C                          contained in the xy plane 
C REAL*8  ARMUNI         : Conversion factors for the charge density
C **********************************************************************

C Internal variables ---------------------------------------------------

      CHARACTER 
     .  OGP*22, OGP_DEFECT*22,
     .  CPF*22, CPF_DEFECT*22,
     .  UCD*22, UCD_DEFECT*22

      INTEGER
     .  IUNIT, IX, JX, NPX_DEFECT, NPY_DEFECT, NPZ_DEFECT,
     .  IND1, IND2, IND3

      DOUBLE PRECISION
     .  ORIGIN(3), XDIR(3)

      LOGICAL 
     .  LEQI, COLIN

      EXTERNAL  LEQI

      DATA ORIGIN /0.D0,0.D0,0.D0/
      DATA XDIR   /1.D0,0.D0,0.D0/
      DATA IND1   /1/
      DATA IND2   /2/
      DATA IND3   /3/
      DATA COLIN  /.FALSE./


      CPF_DEFECT = 'Bohr'

      CPF = FDF_STRING('Denchar.CoorUnits',CPF_DEFECT)

      IF (LEQI(CPF,'bohr')) then
        ISCALE = 1
      ELSEIF (LEQI(CPF,'ang')) then
        ISCALE = 2
      ELSE
        WRITE(6,*)'readpla: ERROR Denchar.CoorUnits must be Ang or Bohr'
        STOP
      ENDIF

      UCD_DEFECT = 'Ele/bohr**3'
      UCD = FDF_STRING('Denchar.DensityUnits',UCD_DEFECT)
      IF (LEQI(UCD,'ele/bohr**3')) then
        IUNITCD = 1
      ELSEIF (LEQI(UCD,'ele/ang**3')) then
        IUNITCD = 2
      ELSEIF (LEQI(UCD,'ele/unitcell')) then
        IUNITCD = 3
      ELSE
       WRITE(6,'(A)')' readpla: ERROR   Wrong Option in Units of      '
       WRITE(6,'(A)')' readpla:  Charge Density                       '
       WRITE(6,'(A)')' readpla:  You must choose one of the following:' 
       WRITE(6,'(A)')' readpla:                                       '
       WRITE(6,'(A)')' readpla:      - Ele/bohr**3                    '
       WRITE(6,'(A)')' readpla:      - Ele/ang**3                     '
       WRITE(6,'(A)')' readpla:      - Ele/unitcell                   '
       STOP
      ENDIF

      NPX_DEFECT = 50
      NPY_DEFECT = 50
      NPZ_DEFECT = 50
      NPX = FDF_INTEGER('Denchar.NumberPointsX',NPX_DEFECT)
      NPY = FDF_INTEGER('Denchar.NumberPointsY',NPY_DEFECT)
      IF (IDIMEN .EQ. 2) THEN
        NPZ = 1
      ELSE IF (IDIMEN .EQ. 3) THEN
        NPZ = FDF_INTEGER('Denchar.NumberPointsZ',NPZ_DEFECT)
      ENDIF

      XMIN = FDF_PHYSICAL('Denchar.MinX',-3.D0,'Bohr')
      XMAX = FDF_PHYSICAL('Denchar.MaxX', 3.D0,'Bohr')
      YMIN = FDF_PHYSICAL('Denchar.MinY',-3.D0,'Bohr')
      YMAX = FDF_PHYSICAL('Denchar.MaxY', 3.D0,'Bohr')
      IF (IDIMEN .EQ. 2) THEN
        ZMIN = 0.D0
        ZMAX = 0.D0
      ELSE IF (IDIMEN .EQ. 3) THEN
        ZMIN = FDF_PHYSICAL('Denchar.MinZ',-3.D0,'Bohr')
        ZMAX = FDF_PHYSICAL('Denchar.MaxZ', 3.D0,'Bohr')
      ENDIF

      OGP_DEFECT = 'NormalVector'
      OGP = FDF_STRING('Denchar.PlaneGeneration',OGP_DEFECT)
      IF (LEQI(OGP,'normalvector')) then
        IOPTION = 1
      ELSEIF (LEQI(OGP,'twolines')) then
        IOPTION = 2
      ELSEIF (LEQI(OGP,'threepoints')) then
        IOPTION = 3
      ELSEIF (LEQI(OGP,'threeatomicindices')) then
        IOPTION = 4
      ELSE
       WRITE(6,'(A)')' readpla: ERROR Wrong Option to Generate Plane  '
       WRITE(6,'(A)')' readpla:  You must choose one of the following:'
       WRITE(6,'(A)')' readpla:                                       '
       WRITE(6,'(A)')' readpla:      - NormalVector                   '
       WRITE(6,'(A)')' readpla:      - TwoLines                       '
       WRITE(6,'(A)')' readpla:      - ThreePoints                    '
       WRITE(6,'(A)')' readpla:      - ThreeAtomicIndices             '
       STOP
      ENDIF

      IF ( FDF_BLOCK('Denchar.CompNormalVector',IUNIT) ) THEN
        READ(IUNIT,*)(NORMAL(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.Comp2Vectors',IUNIT) ) THEN
        READ(IUNIT,*)(DIRVER1(IX),IX=1,3)
        READ(IUNIT,*)(DIRVER2(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.Coor3Points',IUNIT) ) THEN
        READ(IUNIT,*)(COORPO(1,IX),IX=1,3)
        READ(IUNIT,*)(COORPO(2,IX),IX=1,3)
        READ(IUNIT,*)(COORPO(3,IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.Indices3Atoms',IUNIT) ) THEN
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
 

      IF ( FDF_BLOCK('Denchar.PlaneOrigin',IUNIT) ) THEN
        READ(IUNIT,*)(ORIGIN(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.X_Axis',IUNIT) ) THEN
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

      END subroutine readpla

C!-----------------------------------------------------------------

      SUBROUTINE READWAVES(NSPIN,NORB,IFLAG,NWF,PSI,E,IND)

C Reads the wavefunctions and energies from a file written by Siesta
C
C P. Ordejon, July 2003
C **************** INPUT ********************************************
C INTEGER NSPIN     : Number of spin components
C INTEGER NORB      : Number of basis orbitals
C INTEGER IFLAG     : 0=only read and return number of wavefunctions
C                     1=actually read wavefunctions
C **************** INPUT OR OUTPUT **********************************
C INTEGER NWF       : Number of wavefunctions to read 
C                     input(output) if IFLAG=0(1)
C **************** OUTPUT *******************************************
C REAL*8 PSI(NORB,NWF,NSPIN): Wavefunctions
C REAL*8 E(NWF,NSPIN)       : Eigenvalues
C INTEGER IND(NWF)          : List of indexes of wavefunctions
C *******************************************************************

C Modules

      use fdf


      IMPLICIT NONE

      INTEGER IFLAG, NSPIN, NORB, NWF
      INTEGER IND(NWF)
      DOUBLE PRECISION PSI(NORB,NWF,NSPIN), E(NWF,NSPIN)

C INTERNAL VARIABLES .............
      INTEGER UNIT, NK, NSP, NUO, ISPIN, IISPIN, IWF, IIWF, IORB
      INTEGER IDUMB
      DOUBLE PRECISION REPSI,IMPSI

      CHARACTER PASTE*33
      CHARACTER, SAVE :: SNAME*30, FNAME*33
      CHARACTER CHDUMB*20

      SAVE UNIT
      EXTERNAL PASTE
C ..................


c      write(6,*) NORB,NWF,NSPIN

c      write(6,*) 'In readwaves with iflag=',iflag

      IF (IFLAG .EQ. 0) THEN
        SNAME = FDF_STRING('SystemLabel','siesta')
        FNAME = PASTE(SNAME,'.WFS')

        CALL IO_ASSIGN(UNIT)
        OPEN (UNIT, FILE=FNAME, FORM='unformatted', STATUS='unknown')

c        write(6,*) 'opening unit=', unit

        READ(UNIT) NK
        IF (NK .NE. 1) THEN
          WRITE(6,*) 'Wavefunctions file contains more then 1 k-point'
          WRITE(6,*) 'DENCHAR can only handle the Gamma point!!'
          STOP
        ENDIF
        READ(UNIT) NSP
        IF (NSP .NE. NSPIN) THEN
          WRITE(6,*) 'NSPIN is not consistent between data files!'
          STOP
        ENDIF
        READ(UNIT) NUO
        IF (NUO .NE. NORB) THEN
          WRITE(6,*) 'Nr. of orbs is not consistent between data files!'
          STOP
        ENDIF

        READ(UNIT) 
        READ(UNIT) IISPIN
        READ(UNIT)NWF

        REWIND(UNIT)
        
c        write(6,*) 'Exiting iflag= ',iflag
         
        
        RETURN

      ELSE IF (IFLAG .EQ. 1) THEN

c        write(6,*) 'trying to read in unit=', unit
        READ(UNIT) NK
        IF (NK .NE. 1) THEN
          WRITE(6,*) 'Wavefunctions file contains more then 1 k-point'
          WRITE(6,*) 'DENCHAR can only handle the Gamma point!!'
          STOP
        ENDIF
        READ(UNIT) NSP
        IF (NSP .NE. NSPIN) THEN
          WRITE(6,*) 'NSPIN is not consistent between data files!'
          STOP
        ENDIF
        READ(UNIT) NUO
        IF (NUO .NE. NORB) THEN
          WRITE(6,*) 'Nr. of orbs is not consistent between data files!'
          STOP
        ENDIF

        DO ISPIN = 1,NSPIN
          READ(UNIT) 
          READ(UNIT) IISPIN
c          write(6,*) 'spin=',iispin
          IF (IISPIN .NE. ISPIN) THEN
            WRITE(6,*) 'Inconsistent order of spins in wavefuncs. file!'
            STOP
          ENDIF
          READ(UNIT)NWF
c          write(6,*) 'nwfs=',nwf
c          write(6,*) 'number of orbitals=',norb
 

          DO IWF=1,NWF
            READ(UNIT) IND(IWF)
c            write(6,*) 'orbital index ',iwf,' = ',ind(iwf)
            READ(UNIT) E(IWF,ISPIN)
c            write(6,*) 'energy = ',e(iwf,ispin)
            DO IORB = 1, NORB
              READ(UNIT) IDUMB,CHDUMB,IDUMB,IDUMB,IDUMB,CHDUMB,
     .                   REPSI,IMPSI
c              write(6,*) repsi,impsi
              IF (DABS(IMPSI) .GT. 1.0D-10) 
     .          WRITE(6,*) 'Warning: complex wavefunctions in file!'
              PSI(IORB,IWF,ISPIN)=REPSI
            ENDDO
          ENDDO
        ENDDO

        CLOSE (UNIT)
        CALL IO_CLOSE(UNIT)
c        write(6,*) 'Exiting iflag= ',iflag
      ELSE
        WRITE(6,*) 'IFLAG must be either 0 or 1 in READWAVE!!'
        STOP
      ENDIF
           

      END subroutine readwaves
        
C!-----------------------------------------------------------------------

      SUBROUTINE WROUT(IDIMEN, CHARGE, WAVES, IOPTION, NORMAL, COORPO, 
     .                 DIRVER1, DIRVER2, 
     .                 NPX, NPY, NPZ, XMIN, XMAX, YMIN, YMAX, 
     .                 ZMIN, ZMAX, IUNITCD,
     .                 MAXATOM, NAPLA, INDICES, XAPLA )

C **********************************************************************
C Dump input data to ouput
C Modified to make general writeout is done only to standard output, 
C not to the individual data files.
C
C Written by J. Junquera Feb '99
C Modified by P. Ordejon, June 2003
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      LOGICAL
     .  CHARGE, WAVES

      INTEGER
     .  IOPTION,  NPX, NPY, NPZ, IUNITCD, MAXATOM, NAPLA, IDIMEN,
     .  INDICES(MAXATOM)

      DOUBLE PRECISION
     .  NORMAL(3), COORPO(3,3), DIRVER1(3), DIRVER2(3), XAPLA(3,MAXATOM)

      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX

C **************  INPUT  ***********************************************
C INTEGER IDIMEN         : 2D or 3D run
C LOGICAL CHARGE         : Are we writting charge output?
C LOGICAL WAVES          : Are we writting wavefunctions output?
C INTEGER IDIMEN         : 2D or 3D run
C INTEGER IOPTION        : Option to generate the plane
C                          1 = Normal vector
C                          2 = Two vectors belonging to the plane

C                          3 = Three points of the plane
C                          4 = Three atomic indices
C REAL*8  NORMAL(3)      : Components of the normal vector
C REAL*8  COORPO(3,3)    : Coordinates of the three points used to define
C                          the plane 
C REAL*8  DIRVER(3)      : Components of two vector contained in the plane
C INTEGER NPX,NPY,NPZ    : Number of points along x and y and z 
C REAL*8  XMIN, XMAX     : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX     : Limits of the plane in the PRF for y-direction
C REAL*8  ZMIN, ZMAX     : Limits of the grid z-direction
C INTEGER IUNITCD        : Unit of the charge density
C INTEGER MAXATOM        : Total number of atoms in supercell
C INTEGER NAPLA          : Number of atoms whose coordiantes has been rotated   
C INTEGER INDICES(MAXATOM): Indices of tha atoms whose coordinates has 
C                           been roated
C REAL*8  XAPLA(3,MAXATOM): Atomic coordiantes in the in-plane reference frame
C **********************************************************************

C ***************  INTERNAL VARIABLES **********************************
      CHARACTER
     .  SNAME*30

      INTEGER
     .  IX, IP, IA, UNIT1

      LOGICAL, SAVE :: FRSTME

      DATA FRSTME /.TRUE./



C Open files to store charge density -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')

      UNIT1 = 6


C Write general information only if called for the first time

      IF (FRSTME) THEN


        WRITE(UNIT1,'(A)')
     .    '                          ************************       '
        WRITE(UNIT1,'(A)')
     .    '                          *  WELCOME TO DENCHAR  *       '
        WRITE(UNIT1,'(A)')
     .    '                          ************************       '

        WRITE(UNIT1,'(A,A)')
     .    '  You are running DENCHAR for system: ',SNAME
        WRITE(UNIT1,'(A)')
     .    '  '

        WRITE(UNIT1,'(A)')
        IF (IDIMEN .EQ. 2) THEN
          WRITE(UNIT1,'(A)')
     .   '  You have chosen the 2D mode. Values of the functions'
          WRITE(UNIT1,'(A)')
     .    '  will be given in a 2D grid'
        ELSE IF (IDIMEN .EQ. 3) THEN
          WRITE(UNIT1,'(A)')
     .   '  You have chosen the 3D mode. Values of the functions'
          WRITE(UNIT1,'(A)')
     .    '  will be given in a 3D grid, in Gaussian Cube format'
        ENDIF

        WRITE(UNIT1,'(A)')
        WRITE(UNIT1,'(A,/,A,I5)')
     .    '  Number of points in the x-direction : ',
     .    '  ', NPX
        WRITE(UNIT1,'(A,/,A,I5)')
     .    '  Number of points in the y-direction : ',
     .    '  ', NPY
        IF (IDIMEN .EQ. 3)
     .    WRITE(UNIT1,'(A,/,A,I5)')
     .    '  Number of points in the z-direction : ',
     .    '  ', NPZ
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Minimum value of the x-component of the window : ',
     .    '  ', XMIN,' bohrs'
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Maximum value of the x-component of the window : ',
     .    '  ', XMAX,' bohrs'
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Minimum value of the y-component of the window : ',
     .    '  ', YMIN,' bohrs'
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Maximum value of the y-component of the window : ',
     .    '  ', YMAX,' bohrs'
        IF (IDIMEN .EQ. 3) THEN
          WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Minimum value of the z-component of the window : ',
     .    '  ', ZMIN,' bohrs'
          WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Maximum value of the z-component of the window : ',
     .    '  ', ZMAX,' bohrs'
        ENDIF

        WRITE(UNIT1,'(A)')
     .    '  '
        WRITE(UNIT1,'(A,/A)')
     .    '  The options you have chosen to generate the plane',
     .    '  are the following: '

        IF( IOPTION .EQ. 1 ) THEN

          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : NormalVector'
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Components of the normal vector : ',
     .    '  ',(NORMAL(IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Origin of the plane : ',
     .    '  ',(COORPO(1,IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Another point to define the X direction : ',
     .    '  ',(COORPO(2,IX),IX=1,3)

        ELSEIF( IOPTION .EQ. 2 ) THEN 
        
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : TwoLines'
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Components of the first vector inside the plane :',
     .    '  ',(DIRVER1(IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Components of the second vector inside the plane:',
     .    '  ',(DIRVER2(IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Origin of the plane : ',
     .    '  ',(COORPO(1,IX),IX=1,3)

        ELSEIF( IOPTION .EQ. 3 ) THEN 

          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : ThreePoints'
          WRITE(UNIT1,'(A)')
     .    '  Coordinates of three points in the plane : '
          DO IP = 1,3
            WRITE(UNIT1,'(A,3F12.5)')
     .    '  '  ,(COORPO(IP,IX),IX=1,3)
          ENDDO
      
        ELSEIF( IOPTION .EQ. 4 ) THEN 

          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : ThreeAtomicIndices'
          WRITE(UNIT1,'(A)')
     .    '  Position of the three atoms : '
          DO IP = 1,3
            WRITE(UNIT1,'(A,3F12.5)')
     .    '  '  ,(COORPO(IP,IX),IX=1,3)
          ENDDO

        ENDIF

        IF ( IUNITCD .EQ. 1) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A,/,A)')
     .    '  Unit of the charge density in output files : ',
     .    '  Electrons/(bohr**3)'
        ELSEIF ( IUNITCD .EQ. 2) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A,/,A)')
     .    '  Unit of the charge density in output files : ',
     .    '  Electrons/(angstrom**3)'
        ELSEIF( IUNITCD .EQ. 3) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A,/,A)')
     .    '  Unit of the charge density in output files : ',
     .    '  Electrons/unit cell'
        ENDIF


        IF( NAPLA .NE. 0) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Atomic coordinates in the in-plane reference frame'
          WRITE(UNIT1,'(A,19(1H ),A)')
     .    '  Atomic Index','Atomic coordinates'
          DO IA = 1, NAPLA
            WRITE(UNIT1,'(A,I14,5X,3F15.4)')
     .      ' ',INDICES(IA), (XAPLA(IX,INDICES(IA)),IX=1,3)
          ENDDO
        ENDIF

        FRSTME = .FALSE.

      ENDIF

      IF (CHARGE) THEN
        
        WRITE(UNIT1,'(A)')
     .    '  '
        WRITE(UNIT1,'(A)')
     .    '  You are now computing charge density on the grid'

      ENDIF

      IF (WAVES) THEN
        
        WRITE(UNIT1,'(A)')
     .    '  '
        WRITE(UNIT1,'(A)')
     .    '  You are now computing Wave Functions on the grid'

      ENDIF


      END subroutine wrout

      end module m_denchar_io

      



