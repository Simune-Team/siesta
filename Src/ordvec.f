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
      EXTERNAL          CHKDIM, IORDER, ORDER, ORDIX

      INTEGER MAXAUX
      PARAMETER ( MAXAUX = 1000 )
      INTEGER IAUX(MAXAUX), IV, IV0, IV1, IX, JX

      DO 20 IV = 1,NV
        INDEX(IV) = IV
  20  CONTINUE
      DO 80 IX = NX,1,-1
        IV0 = 0
  30    CONTINUE
          DO 50 IV1 = IV0+1,NV-1
            DO 40 JX = IX+1,NX
              IF (ABS(V(JX,IV1+1)-V(JX,IV1)) .GT. TOL) GOTO 60
  40        CONTINUE
  50      CONTINUE
            IV1 = NV
  60      CONTINUE
          IF (IV1 .GT. IV0+1) THEN
            CALL CHKDIM( 'ORDVEC', 'MAXAUX', MAXAUX, IV1-IV0, 1 )
            CALL ORDIX(  V(IX,IV0+1), NX, IV1-IV0, IAUX )
            CALL ORDER(  V(1,IV0+1),  NX, IV1-IV0, IAUX )
            CALL IORDER( INDEX(IV0+1), 1, IV1-IV0, IAUX )
          ENDIF
          IV0 = IV1
        IF (IV0 .LT. NV-1) GO TO 30
  80  CONTINUE
      END

