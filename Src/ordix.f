      SUBROUTINE ORDIX( X, M, N, INDX )
C *******************************************************************
C Makes an index table of array X, with size N and stride M.
C Ref: W.H.Press et al. Numerical Recipes, Cambridge Univ. Press.
C Adapted by J.M.Soler from routine INDEXX of Num. Rec. May'96.
C *************** INPUT *********************************************
C REAL*8  X(M,N)   : Array with the values to be ordered
C INTEGER M, N     : Dimensions of array X
C *************** OUTPUT ********************************************
C INTEGER INDEX(N) : Array which gives the increasing order of X(1,I):
C                    X(1,INDEX(I)) .LE. X(1,INDEX(I+1)) )
C *************** USAGE *********************************************
C Example to order atomic positions X(I,IA), I=1,3, IA=1,NA by
C increasing z coordinate:
C    CALL ORDIX( X(3,1), 3, NA, INDEX )
C    CALL ORDER( X(1,1), 3, NA, INDEX )
C *******************************************************************
      IMPLICIT          NONE
      INTEGER           I, N, INDX(N), INDXT, IR, J, L, M
      DOUBLE PRECISION  X(M,N), Q

      DO 1 J=1,N
         INDX(J)=J
   1  CONTINUE
      IF (N.LE.1) RETURN
      L=N/2+1
      IR=N
   2  CONTINUE
         IF (L.GT.1) THEN
            L=L-1
            INDXT=INDX(L)
            Q=X(1,INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=X(1,INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF (IR.EQ.1) THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
   3     IF (J.LE.IR) THEN
            IF (J.LT.IR) THEN
               IF (X(1,INDX(J)).LT.X(1,INDX(J+1))) J=J+1
            ENDIF
            IF (Q.LT.X(1,INDX(J))) THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GO TO 3
         ENDIF
         INDX(I)=INDXT
      GO TO 2
      END
