C $Id: superx.f,v 1.3 1999/02/21 00:17:47 emilio Exp $

      SUBROUTINE SUPERX( UCELL, NSC, NA, MAXA, XA, SCELL )

C **********************************************************************
C Generates supercell vectors and atomic positions.
C Written by J.M.Soler, August 1998
C *************** Input ************************************************
C Real*8  UCELL(3,3)  : Unit cell vectors UCELL(Ixyz,Ivector)
C Integer NSC(3)      : Number of cells in each supercell direction:
C                         SCELL(ix,i) = UCELL(ix,i) * NSC(i)
C Integer NA          : Number of atoms in unit cell
C Integer MAXA        : Second dimension of XA
C *************** Input and output *************************************
C Real*8  XA(3,MAXA)  : Atomic positions in unit cell (input) and
C                       in supercell (output), in cartesian coord.
C Real*8  SCELL(3,3)  : Supercell vectors
C *********** Units ****************************************************
C Units of CELL and XA are arbitrary but must be the same
C *********** Behavior *************************************************
C - If NA*NCELLS > MAXA (where NCELLS is the total number of cells),
C   the supercell atomic coordinates are not generated.
C - The first supercell atoms are those of the initial unit cell, i.e.
C   the positions XA(i,ia) for (ia.le.NA) are not modified.
C - The remaining atoms are ordered by unit cells, i.e. the atom ia
C   is equivalent to the unit-cell atom ja=MOD(ia-1,NA)+1
C **********************************************************************

      IMPLICIT          NONE
      INTEGER           MAXA, NA, NSC(3)
      DOUBLE PRECISION  SCELL(3,3), UCELL(3,3), XA(3,MAXA)

C Internal variables
      INTEGER           I, I1, I2, I3, IA, IX, JA, NCELLS
      DOUBLE PRECISION  XC(3)

C Find supercell vectors
      DO 10 I = 1,3
        DO 5 IX = 1,3
          SCELL(IX,I) = UCELL(IX,I) * NSC(I)
    5   CONTINUE
   10 CONTINUE

C Expand atomic positions to supercell
      NCELLS = NSC(1) * NSC(2) * NSC(3)
      IF (NA*NCELLS .LE. MAXA) THEN
        IA = 0
        DO 60 I3 = 0,NSC(3)-1
        DO 50 I2 = 0,NSC(2)-1
        DO 40 I1 = 0,NSC(1)-1
          DO 15 IX = 1,3
            XC(IX) = UCELL(IX,1)*I1 + UCELL(IX,2)*I2 + UCELL(IX,3)*I3
   15     CONTINUE
          DO 30 JA = 1,NA
            IA = IA + 1
            DO 20 IX = 1,3
              XA(IX,IA) = XA(IX,JA) + XC(IX)
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
   50   CONTINUE
   60   CONTINUE
      ENDIF

      END
