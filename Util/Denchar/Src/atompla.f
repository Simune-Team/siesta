      SUBROUTINE ATOMPLA( NA, ORIGIN, XA, MROT, NATINPLA, INDICES, 
     .                    ISCALE, XAPLANE )
C **********************************************************************
C Projection of the coordinates of some selected atoms on the 
C plane-reference-frame. Atoms in plane will have the third coordinate 
C equal to zero.
C Coded by J. Junquera May '99
C **********************************************************************

      IMPLICIT NONE

      INCLUDE 'fdfdefs.h'

      INTEGER
     .  NA, NATINPLA, INDICES(NA), ISCALE

      DOUBLE PRECISION
     .  ORIGIN(3), XA(3,NA), MROT(3,3), XAPLANE(3,NA)

C ******* INPUT ********************************************************
C INTEGER NA             : Number of atoms
C REAL*8  ORIGIN(3)      : Origin of the plane reference frame
C REAL*8  XPLANE(3)      : X-vector of the in-plane reference frame 
C REAL*8  YPLANE(3)      : Y-vector of the in-plane reference frame
C REAL*8  NORMAL(3)      : Normal vector of the plane 
C REAL*8  XA(3,NA)       : Atomic coordinates in lattice reference frame
C REAL*8  MROT(3,3)      : Rotation matrix from the lattice reference frame 
C                          to the in-plane reference frame
C INTEGER ISCALE         : Units of the atomic coordinates
C                          ISCALE = 1 => bohrs, ISCALE = 2 => Ang
C ******* OUTPUT *******************************************************
C INTEGER NATINPLA       : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA)    : Atomic indices of the atoms whose coordinates
C                          will be rotated
C REAL*8  XAPLANE(3,NA)  : Atomic coordinates in plane reference frame
C **********************************************************************

      CHARACTER
     .  LINE*130, NAMES*80

      INTEGER 
     .  NATINPL_DEFECT, IUNIT, IAT, LASTC, NN, NV, NI,
     .  LC(0:3), INTEGS(4), NR, IX
 
      DOUBLE PRECISION
     .  VALUES(4), REALS(4), VAUX1(3), VAUX2(3)

      LOGICAL 
     .  ATINPLA

      EXTERNAL 
     .  MATVECT, PARSE 

C **********************************************************************
C INTEGER NATINPLA       : Number of atoms whose coordinates will be 
C                          rotated from the lattice reference frame to 
C                          the in-plane reference frame
C INTEGER NN             : Number of names present in line
C INTEGER LC(0:3)        : Last character of each name in names string.
C                          Notice that first index is 0.
C CHARACTER NAMES*80     : Non-number strings present in line
C INTEGER NV             : Number of values present in line
C REAL*8 VALUES(4)       : Numbers present in line(integer or real)
C INTEGER NI             : Number of integers present in line
C INTEGER INTEGS(4)      : Integer numbers present in line
C INTEGER NR             : Number of reals present in line
C REAL*8  REALS(4)       : Real numbers present in line
C REAL*8  VAUX1(3)       : Auxiliar vector
C REAL*8  VAUX2(3)       : Auxiliar vector
C LOGICAL ATINPLA        : Is the block 2D.AtomsInPlane present in fdf?
C **********************************************************************

C Read fdf data block '2D.AtomsInPlane' --------------------------------
      NATINPL_DEFECT = 0

      ATINPLA = FDF_BLOCK('2D.AtomsInPlane',IUNIT)

      IF ( ATINPLA ) THEN

        NATINPLA = 0
        DO IAT = 1, NA + 1
          READ(IUNIT,'(A)',END=50) LINE
          LASTC = INDEX(LINE,'#') - 1
          IF (LASTC .LE. 0) LASTC = LEN(LINE)
          CALL PARSE( LINE(1:LASTC), NN, LC, NAMES, NV, VALUES,
     .                NI, INTEGS, NR, REALS )    
          IF (NN .GE. 1 .AND. NAMES(LC(0)+1:LC(1)).EQ.'%endblock') THEN
C         End data reading 
          GOTO 50
          ELSEIF (NI .EQ. 1) THEN
             NATINPLA = NATINPLA + 1
             INDICES(NATINPLA) = INTEGS(1)
          ELSE
C            Print bad sintax error and stop
             GOTO 40
          ENDIF
        ENDDO

        WRITE(6,'(A)')
     .    'ATOMPLA: ERROR: Too many atom entries in 2D.AtomsInPlane'
        STOP 'ATOMPLA: ERROR: Too many atom entries in 2D.AtomsInPlane'

 40     CONTINUE
        WRITE(6,*)
     .    'ATOMPLA: ERROR: BAD SYNTAX IN 2D.AtomsInPlane, line ', IAT
        STOP 'ATOMPLA: ERROR: BAD SYNTAX IN 2D.AtomsInPlane'

 50     CONTINUE
      ENDIF

C Rotate the coordinates -----------------------------------------------
        DO IAT = 1, NATINPLA

          DO IX = 1,3
            VAUX1(IX) = XA(IX,INDICES(IAT)) - ORIGIN(IX)
          ENDDO

          CALL MATVECT(MROT, VAUX1, VAUX2)

          DO IX = 1,3
            XAPLANE(IX,INDICES(IAT)) = VAUX2(IX)
          ENDDO
  
        IF (ISCALE .EQ. 2) THEN
          DO IX = 1,3
           write(6,*)xaplane(ix,indices(iat))
           XAPLANE(IX,INDICES(IAT))=XAPLANE(IX,INDICES(IAT))*0.529177D0       
           write(6,*)xaplane(ix,indices(iat))
          ENDDO
        ENDIF

        ENDDO

        RETURN
        END
