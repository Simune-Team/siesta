      module sorting

      CONTAINS

      SUBROUTINE ORDVEC( TOL, NX, NV, V, INDEX )
      
C **********************************************************************
C Orders a set of vectors in a well defined order: by last coordinate;
C if last coordinate is equal, by before-last corrdinate, etc.
C Written by J.Soler. October 1997.
C ************ INPUT ***************************************************
C REAL*8  TOL      : Tolerance to consider two coordinates equal.
C INTEGER NX       : Number of vector coordinates
C INTEGER NV       : Number of vectors
C REAL*8  V(NX,NV) : Vectors to be ordered
C ************ OUTPUT **************************************************
C REAL*8  V(NX,NV) : Ordered vectors
C INTEGER INDEX(NV): Order index such that Vout(IX,IV)=Vin(IX,INDEX(IV))
C **********************************************************************

      IMPLICIT          NONE
      INTEGER           NV, NX
      INTEGER           INDEX(NV)
      DOUBLE PRECISION  TOL, V(NX,NV)

      INTEGER           IV, IV0, IV1, IX, JX

      INTEGER :: IAUX(NV)           ! automatic array

      DO IV = 1,NV
        INDEX(IV) = IV
      ENDDO
      DO IX = NX,1,-1
        IV0 = 0
  30    CONTINUE
          DO IV1 = IV0+1,NV-1
            DO JX = IX+1,NX
              IF (ABS(V(JX,IV1+1)-V(JX,IV1)) .GT. TOL) GOTO 60
            ENDDO
          ENDDO
          IV1 = NV
  60      CONTINUE
          IF (IV1 .GT. IV0+1) THEN
            CALL ORDIX(  V(IX:IX,IV0+1:), NX, IV1-IV0, IAUX )
            CALL ORDER(  V(1,IV0+1),  NX, IV1-IV0, IAUX )
            CALL IORDER( INDEX(IV0+1), 1, IV1-IV0, IAUX )
          ENDIF
          IV0 = IV1
        IF (IV0 .LT. NV-1) GO TO 30
      ENDDO

      END subroutine ordvec

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
      INTEGER           I, N, INDX(:), INDXT, IR, J, L, M
      DOUBLE PRECISION  X(:,:), Q

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

      END subroutine ordix


      SUBROUTINE ORDER( X, M, N, INDEX )
C *******************************************************************
C Orders array X(M,N) according to array INDEX(N), which may be
C   generated by routine ORDIX.
C Written by J.M.Soler. May'96.
C *************** INPUT *********************************************
C INTEGER M, N     : Dimensions of array X
C INTEGER INDEX(N) : Array which gives the desired order
C *************** INPUT AND OUTPUT **********************************
C REAL*8  X(M,N) : Array(s) to be ordered: Xout(I,J) = Xin(I,INDEX(J))
C *******************************************************************
      IMPLICIT          NONE
      INTEGER           I, N, INDEX(N), IORDER, ISTORE, J, M
      DOUBLE PRECISION  X(M,N), XI

      DO 40 J = 1,M
        DO 20 I = 1,N
          XI = X(J,I)
          IORDER = I
   10     CONTINUE
          ISTORE = INDEX(IORDER)
          IF (ISTORE .GT. 0) THEN
            IF (ISTORE .EQ. I) THEN
              X(J,IORDER) = XI
            ELSE
              X(J,IORDER) = X(J,ISTORE)
            ENDIF
            INDEX(IORDER) = -INDEX(IORDER)
            IORDER = ISTORE
            GOTO 10
          ENDIF
   20   CONTINUE
        DO 30 I = 1,N
          INDEX(I) = -INDEX(I)
   30   CONTINUE
   40 CONTINUE
      END subroutine order


      SUBROUTINE IORDER( IA, M, N, INDEX )
C *******************************************************************
C Orders integer array IA(M,N) according to array INDEX(N), which
C   may be generated by routine ORDIX.
C Written by J.M.Soler. May'96 and Oct'97.
C *************** INPUT *********************************************
C INTEGER M, N     : Dimensions of array IA
C INTEGER INDEX(N) : Array which gives the desired order
C *************** INPUT AND OUTPUT **********************************
C REAL*8  IA(M,N): Array(s) to be ordered: IAout(I,J) = IAin(I,INDEX(J))
C *******************************************************************
      IMPLICIT NONE
      INTEGER  M, N
      INTEGER  IA(M,N), INDEX(N)

      INTEGER  I, IAI, IORD, ISTORE, J

      DO 40 J = 1,M
        DO 20 I = 1,N
          IAI = IA(J,I)
          IORD = I
   10     CONTINUE
          ISTORE = INDEX(IORD)
          IF (ISTORE .GT. 0) THEN
            IF (ISTORE .EQ. I) THEN
              IA(J,IORD) = IAI
            ELSE
              IA(J,IORD) = IA(J,ISTORE)
            ENDIF
            INDEX(IORD) = -INDEX(IORD)
            IORD = ISTORE
            GOTO 10
          ENDIF
   20   CONTINUE
        DO 30 I = 1,N
          INDEX(I) = -INDEX(I)
   30   CONTINUE
   40 CONTINUE
      END subroutine iorder

      end module sorting
