C $Id: listsc.f,v 1.3 1999/01/31 11:20:06 emilio Exp $

      SUBROUTINE LISTSC( NSC, NO, MAXO, MAXNO, NUMNO, LISTNO )

C *********************************************************************
C Expands neighbour lists from unit cell to supercell
C Written by J.M.Soler, August 1998
C ************** INPUT ************************************************
C INTEGER NSC(3) : Number of cells in each supercell direction:
C                    SuperCell(ix,i) = CELL(ix,i) * NSC(i)
C INTEGER NO     : Number of orbitals in unit cell
C INTEGER MAXO   : Last dimension of NUMNO and LISTNO
C INTEGER MAXNO  : First dimension of LISTNO
C ************** INPUT (IO.LE.NO) AND OUTPUT (IO.GT.NO) ***************
C INTEGER NUMNO(MAXO)    : Number of elements in each column of LISTNO
C INTEGER LISTNO(MAXNO,MAXO): List of neighbors of orbitals in supercell
C *********************************************************************

      IMPLICIT  NONE
      INTEGER   MAXNO, MAXO, NO
      INTEGER   LISTNO(MAXNO,MAXO), NUMNO(*), NSC(3)
      EXTERNAL  IPACK

C Internal variables and arrays
      INTEGER   I, IC(3), ICELL, IL, IO, ION,
     .          JC(3), JCELL, JO, JON, NCELLS

C Find total number of unit cells in supercell
      NCELLS = NSC(1) * NSC(2) * NSC(3)

C Expand NUMNO and LISTNO to supercell
      IF (NCELLS.GT.1 .AND. NO*NCELLS.LE.MAXO) THEN
        DO ICELL = 1,NCELLS
          CALL IPACK( -1, 3, NSC, IC, ICELL )
          DO IO = 1,NO
            JO = IO + NO * (ICELL-1)
            NUMNO(JO) = NUMNO(IO)
            IF (NUMNO(IO) .LE. MAXNO) THEN
              DO IL = 1,NUMNO(IO)
C                 ION is a neighbour orbital of IO
                ION = LISTNO(IL,IO)
C                 JCELL is the unit cell where ION is
                JCELL = (ION-1) / NO + 1
C                 Find the orbital JON which has, relative to JO, the
C                 same position that ION has relative to IO
                CALL IPACK( -1, 3, NSC, JC, JCELL )
                DO I = 1,3
                  JC(I) = JC(I) + IC(I)
                ENDDO
                CALL IPACK( +1, 3, NSC, JC, JCELL )
C                 Now JCELL is the unit cell where JON is
                JON = 1 + MOD(ION-1,NO) + NO * (JCELL-1)
                LISTNO(IL,JO) = JON
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      END




      SUBROUTINE IPACK( TASK, ND, N, I, IND )

C *********************************************************************
C Packs/unpacks several integer indexes into/out of one
C Written by J.M.Soler. July 1998
C ******** INPUT ******************************************************
C INTEGER TASK : Task switch: TASK=+1 => Pack I into IND
C                             TASK=-1 => Unpack I out of IND
C INTEGER ND   : Number of indexes to pack/unpack (dimension of I)
C INTEGER N(ND): Upper range of I. Must be N(j)>0 for all j.
C ******** INPUT or OUTPUT (Depending of TASK) ************************
C INTEGER I(ND): Indexes to pack/unpack
C INTEGER IND  : Combined (packed) index such that
C                IND = 1 + I(1) + N(1)*I(2) + N(1)*N(2)*I(3) + ...
C                where 0 <= I(j) < N(j)
C ******** BEHAVIOR ***************************************************
C Notice that the range of I(j) begins at 0, and that of IND at 1
C A modulus operation is done to bring I(j) to the range (0:N(j)-1)
C Does not check that N(j) > 0
C If TASK=-1 and IND is not in the range (1:N(1)*N(2)*...*N(ND)),
C   the program stops with an error message
C *********************************************************************

      IMPLICIT NONE
      INTEGER  IND, ND, TASK
      INTEGER  I(ND), N(ND)

      INTEGER  IJ, J

      IF (TASK .GT. 0) THEN
        IND = 0
        DO J = ND,1,-1
          IJ = MOD( I(J)+1000*N(J), N(J) )
          IND = IJ + N(J) * IND
        ENDDO
        IND = 1 + IND
      ELSE
        IF (IND.LT.1) STOP 'IPACK: ERROR: IND < 1'
        IJ = IND - 1
        DO J = 1,ND
          I(J) = MOD( IJ, N(J) )
          IJ = (IJ-I(J)) / N(J)
        ENDDO
        IF (IJ.GT.0) STOP 'IPACK: ERROR: IND out of range'
      ENDIF

      END

