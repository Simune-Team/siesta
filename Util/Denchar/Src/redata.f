
      SUBROUTINE REDATA( NSPIN, NUO, NO, NA, 
     .                   MAXO, MAXA, MAXNO, MAXSPIN, 
     .                   CELL, RMAXO, XA, 
     .                   LASTO, ISA, IPHORB, DATM,
     .                   NUMH, LISTH, INDXUO )
C **********************************************************************
C Read tables calculated in siesta
C Coded by J. Junquera 11/98
C **********************************************************************

      IMPLICIT NONE

      INCLUDE 'atom.h'
      INCLUDE 'fdfdefs.h'

      INTEGER 
     .  NSPIN, MAXO, MAXA, MAXSPIN, MAXNO,
     .  NUO, NO, NA, UNIT1, LASTO(0:MAXA), 
     .  ISA(MAXA), IPHORB(MAXO),
     .  NUMH(MAXO), LISTH(MAXNO,MAXO), INDXUO(MAXO),
     .  ISMAX, NOMAX(NSMAX), NKBMAX(NSMAX),
     .  IB, IL, IS, IA, IO, J, NPOLORBSAVE(LMAXD,NSMAX),
     .  LMXOSAVE(NSMAX), LMXKBSAVE(NSMAX), NZETASAVE(0:LMAXD,NSMAX)


      DOUBLE PRECISION
     .  CELL(3,3), XA(3,MAXA), RMAXO, DATM(MAXO)

      DOUBLE PRECISION
     .  TABLE((NTBMAX+2),-(LMAXD+1):NZETMX*(LMAXD+1),NSMAX),
     .  TAB2(NTBMAX,-(LMAXD+1):NZETMX*(LMAXD+1),NSMAX),
     .  TABPOL((NTBMAX+2),LMAXD*NZETMX,NSMAX),
     .  TAB2POL(NTBMAX,LMAXD*NZETMX,NSMAX)


      CHARACTER
     .  SNAME*30, FNAME*33, PASTE*33

      EXTERNAL 
     .  IO_ASSIGN, IO_CLOSE, PASTE

      COMMON/CMTAB/TABLE,TABPOL
      COMMON/CMSPLINE/TAB2,TAB2POL
      COMMON/CMZETA/NZETASAVE
      COMMON/CONTROL/ISMAX,NOMAX,NKBMAX
      COMMON/CMPOLORB/NPOLORBSAVE
      COMMON/CMLMXO/LMXOSAVE
      COMMON/CMLMXKB/LMXKBSAVE

C **********************************************************************
C INTEGER MAXO           : Max. total number of basis orbitals
C INTEGER MAXA           : Max. number of atoms
C INTEGER MAXNO          : Max. number of neighbour orbitals
C INTEGER MAXSPIN        : Max. number of spines (1 or 2)
C INTEGER NSPIN          : Number of different spin polarizations
C                          Nspin = 1 => Non polarized, Nspin = 2 => Polarized 
C INTEGER NUO            : Number of orbitals in unit cell
C INTEGER NO             : Total number of atomic orbitals
C INTEGER NA             : Number of atoms 
C REAL*8  CELL(3,3)      : Unit cell vectors CELL(IXYZ,IVECT)
C REAL*8  XA(3,MAXA)     : Atomic position in cartesian coordinates
C REAL*8  RMAXO          : Maximum cutoff radius for atomic orbitals
C INTEGER LASTO(0:MAXA)  : Position of last orbital of each atom
C INTEGER ISA(MAXA)      : Species index of each atom
C INTEGER IPHORB(MAXO)   : Orbital index of each atom in its orbital
C REAL*8  DATM(MAXO)     : Neutral atom charge of each orbital
C INTEGER NUMH(MAXO)     : Control vector of Density Matrix
C                          (number of nonzero elements of each row)
C INTEGER LISTH(MAXNO,MAXO) : Control vector of Density Matrix
C                          (list of nonzero elements of each row)
C INTEGER INDXUO(MAXO)   : Index of equivalent orbital in unit cell
C **********************************************************************

      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = PASTE(SNAME,'.PLD')

      CALL IO_ASSIGN(UNIT1)

      OPEN (UNIT = UNIT1, FILE = FNAME, FORM = 'UNFORMATTED',
     .      STATUS = 'UNKNOWN')
        REWIND(UNIT1)
C Read tables from files -----------------------------------------------

CIn the previous version this loop only was executed from
CIL= -LMAXD+1 instead of -(LMAXD+1).
C I belive that now is right! DSP.
      DO 100 IB = 1, NTBMAX+2
        DO 110 IL = -(LMAXD+1), NZETMX*(LMAXD+1)
          DO 120 IS = 1, NSMAX
            READ(UNIT1)TABLE(IB,IL,IS)
 120      CONTINUE
 110    CONTINUE
 100  CONTINUE

      DO 200 IB = 1, NTBMAX+2
        DO 210 IL = 1, NZETMX*LMAXD
          DO 220 IS = 1, NSMAX
            READ(UNIT1)TABPOL(IB,IL,IS)
 220      CONTINUE
 210    CONTINUE
 200  CONTINUE

CIn the previous version this loop only was executed from
CIL= -LMAXD+1 instead of -(LMAXD+1).
C I belive that now is right! DSP.
      DO 300 IB = 1, NTBMAX
        DO 310 IL = -(LMAXD+1), NZETMX*(LMAXD+1)
          DO 320 IS = 1,NSMAX
            READ(UNIT1)TAB2(IB,IL,IS)
 320      ENDDO
 310    ENDDO
 300  ENDDO


      DO 400 IB = 1, NTBMAX
        DO 410 IL = 1, NZETMX*LMAXD
          DO 420 IS = 1,NSMAX
            READ(UNIT1)TAB2POL(IB,IL,IS)
 420      ENDDO
 410    ENDDO
 400  ENDDO
      
      DO 500 IS = 1,NSMAX
         READ(UNIT1)LMXOSAVE(IS)
 500  ENDDO

      DO 600 IS = 1, NSMAX
         READ(UNIT1)LMXKBSAVE(IS)
 600  ENDDO

      DO 700 IL = 0,LMAXD
        DO 800 IS = 1,NSMAX
          READ(UNIT1)NZETASAVE(IL,IS)
 800    ENDDO
 700  ENDDO
       
      DO IS = 1, NSMAX
        DO IL =1 , LMAXD
          READ(UNIT1)NPOLORBSAVE(IL,IS)
        ENDDO 
      ENDDO 

 
      READ(UNIT1)ISMAX

      DO IS = 1, NSMAX
          READ(UNIT1)NOMAX(IS), NKBMAX(IS)
      ENDDO
     

        DO IL = 1, MAXO
          READ(UNIT1)DATM(IL)
        ENDDO

        READ(UNIT1)NUO
        READ(UNIT1)NO
        READ(UNIT1)NA
        READ(UNIT1)NSPIN

        IF(NO .GT. MAXO) THEN
          WRITE(6,'(a)')'Error in redata.f: NO > MAXO'
          WRITE(6,*)'NO = ',NO,' MAXO = ', MAXO
          WRITE(6,'(a)')'Have you copy the correct .h files into '
          WRITE(6,'(a)')'Src and recompile? '
          STOP
        ENDIF

        IF(NA .GT. MAXA) THEN
          WRITE(6,'(a)')'Error in redata.f: NA > MAXA'
          WRITE(6,*)'NA = ',NA,' MAXA = ', MAXA
          WRITE(6,'(a)')'Have you copy the correct .h files into '
          WRITE(6,'(a)')'Src and recompile? '
          STOP
        ENDIF

        IF(NSPIN .GT. MAXSPIN) THEN
          WRITE(6,'(a)')'Error in redata.f: NSPIN > MAXSPIN'
          WRITE(6,*)'NSPIN = ',NSPIN,' MAXSPIN = ', MAXSPIN
          WRITE(6,'(a)')'Have you copy the correct .h files into '
          WRITE(6,'(a)')'Src and recompile?'
          STOP
        ENDIF

        DO IA = 1,3
          READ(UNIT1)(CELL(J,IA),J=1,3)
        ENDDO

        DO IA = 1,MAXA
          READ(UNIT1)(XA(J,IA),J=1,3),ISA(IA)
        ENDDO

        DO IA = 0, MAXA
          READ(UNIT1)LASTO(IA)
        ENDDO

        READ(UNIT1)RMAXO
        
        DO IO = 1,MAXO
          READ(UNIT1)IPHORB(IO),INDXUO(IO)
        ENDDO
 
        DO IO = 1, MAXO
          READ(UNIT1)NUMH(IO)
          DO J = 1,NUMH(IO)
            READ(UNIT1)LISTH(J,IO)
          ENDDO
        ENDDO


      CALL IO_CLOSE(UNIT1)

      RETURN

      END
