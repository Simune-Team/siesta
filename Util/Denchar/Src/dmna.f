
      SUBROUTINE DMNA( MAXUO, MAXO, MAXND, NUMD, LISTD, LISTDPTR,
     .                 DATM, DSCFNA )
C **********************************************************************
C Form Density Matrix for neutral and isolated atoms from occupations 
C of basis orbitals in free atom
C Coded by J. Junquera 11/98
C Modified by J. Junquera 07/01
C **********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: 
     .  MAXUO, MAXO, MAXND,
     .  NUMD(MAXUO), LISTD(MAXND), LISTDPTR(MAXUO)
 
      DOUBLE PRECISION, INTENT(IN) ::
     .  DATM(MAXO)
 
      DOUBLE PRECISION, INTENT(OUT) ::
     .  DSCFNA(MAXND)
 
C **** INPUT ***********************************************************
C INTEGER MAXUO               : Max. number of basis orbitals in unit cell
C INTEGER MAXO                : Max. number of basis orbitals in supercell
C INTEGER MAXND               : Nonzero Hamiltonian-matrix element
C                               column indexes for each matrix row
C                               For parallel execution, listh contains the
C                               elements for rows that involve
C                               any locally stored orbitals. In the case
C                               where parallelisation is over K points then
C                               the full listh matrix is needed on every
C                               Node.
C INTEGER NUMD(MAXUO)         : Number of non-zero elements of each row of 
C                               the Density Matrix 
C INTEGER LISTD(MAXND)        : Non-zero Density-Matrix element column 
C                               indexes for each matrix row
C INTEGER LISTDPTR(MAXUO)     : Pointer to where each row of listh starts - 1
C                               The reason for pointing to the element before
C                               the first one is so that when looping over the
C                               elements of a row there is no need to shift by
C                               minus one.
C REAL*8  DATM(MAXO)          : Neutral atom charge of each orbital
C **** OUTPUT **********************************************************
C REAL*8  DSCFNA(MAXND)       : Neutral Atom Density Matrix
C **********************************************************************

C Some internal variables ----------------------------------------------

      INTEGER
     .  IO, JO, IN, IND

C Initialize Neutral Atom Density Matrix -------------------------------
      DSCFNA(:) = 0.D0

      DO IO = 1, MAXUO
        DO IN = 1, NUMD(IO)
          IND = LISTDPTR(IO) + IN
          JO = LISTD(IND)
          IF (IO .EQ. JO) THEN
              DSCFNA(IND) = DATM(IO)
          ENDIF
        ENDDO
      ENDDO

      RETURN
 
      END
