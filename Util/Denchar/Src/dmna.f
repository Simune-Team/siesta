
      SUBROUTINE DMNA( NO, NSPIN, MAXO, MAXNO, NUMH, LISTH, INDXUO, 
     .                 DATM, DSCFNA )
C **********************************************************************
C Form Density Matrix for neutral and isolated atoms from occupations 
C of basis orbitals in free atom
C Coded by J. Junquera 11/98
C **********************************************************************

      IMPLICIT NONE

      INTEGER ORBMAX
      PARAMETER( ORBMAX = 10000 )

      INTEGER 
     .  NO, NSPIN, MAXO, MAXNO, INDXUO(MAXO),
     .  NUMH(MAXO), LISTH(MAXNO,MAXO),
     .  IO, IUO, JO, IN, ISPIN
 
      DOUBLE PRECISION
     .  DATM(MAXO),
     .  DSCFNA(MAXNO,MAXO,NSPIN)
 
      DOUBLE PRECISION
     .  DI(ORBMAX)

      EXTERNAL 
     .  CHKDIM

C **********************************************************************
C INTEGER NO                : Total number of orbitals
C INTEGER NSPIN             : Spin polarization
C INTEGER MAXO              : Max. total number of basis orbitals
C INTEGER MAXNO             : Max. total number of neighbours orbitals
C INTEGER INDXUO(MAXO)      : Equivalent otbital in unit cell
C INTEGER NUMH(MAXO)        : Number of non-zero elements of each row of 
C                             the Density Matrix 
C INTEGER LISTH(MAXNO,MAXO) : Non-zero Density-Matrix element column 
C                             indexes for each matrix row
C REAL*8  DATM(MAXO)        : Neutral atom charge of each orbital              
C REAL*8  DSCFNA(MAXNO,MAXO,NSPIN) : Neutral Atom Density Matrix
C **********************************************************************

C Check dimensions -----------------------------------------------------
      CALL CHKDIM( 'DMNA', 'ORBMAX' , ORBMAX , MAXO   , 1 )

C Expand Datm from Unit Cell to SuperCell ------------------------------
      DO IO = 1,NO
        IUO = INDXUO(IO)
        DI(IO) = DATM(IUO)
      ENDDO


      DO IO = 1, NO
C Initialize Neutral Atom Density Matrix -------------------------------
        DO IN = 1, NUMH(IO)
          DO ISPIN = 1, NSPIN
            DSCFNA(IN,IO,ISPIN) = 0.D0
          ENDDO
        ENDDO

        DO IN = 1, NUMH(IO)
          JO = LISTH(IN,IO)
          IF (IO .EQ. JO) THEN
            IF (NSPIN .EQ. 1) THEN

C No spin polarization -------------------------------------------------
              DSCFNA(IN,IO,1) = DI(IO)
            ENDIF

C Spin polarized case is not implemented yet ---------------------------
          ENDIF
        ENDDO
      ENDDO

      RETURN
 
      END
