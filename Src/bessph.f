      DOUBLE PRECISION FUNCTION BESSPH (L,X)
*
*  RETURNS THE SPHERICAL BESSEL FUNCTION JL(X).
*  REF: ABRAMOWITZ AND STEGUN, FORMULAS 10.1.2 AND 10.1.19
*  WRITTEN BY J.SOLER (JSOLER AT EMDUAM11). NOV/89.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TINY=1.D-15,NTERMS=100)
      SWITCH=MAX(1,2*L-1)
      IF (ABS(X).LT.SWITCH) THEN
*        USE POWER SERIES
         TERM=ONE
         DO 10 I=1,L
            TERM=TERM*X/(2*I+1)
   10    CONTINUE
         X2=X*X
         SUM=ZERO
         DO 20 I=1,NTERMS
            SUM=SUM+TERM
            TERM=(-TERM)*X2/(2*I*(2*I+2*L+1))
            IF (ABS(TERM).LT.TINY) GO TO 30
   20    CONTINUE
            WRITE(6,*) 'BESSPH: SERIES HAS NOT CONVERGED. L,X=',L,X
            STOP 3
   30    BESSPH=SUM
      ELSE
*        USE EXPLICIT EXPRESSIONS OR RECURRENCE RELATION
         IF (L.EQ.0) THEN
            BESSPH=SIN(X)/X
         ELSEIF (L.EQ.1) THEN
            BESSPH=(SIN(X)/X-COS(X))/X
         ELSE
            Y=ONE/X
            FNM1=SIN(X)*Y
            FN=(FNM1-COS(X))*Y
            DO 40 N=1,L-1
               FNP1=(2*N+1)*Y*FN-FNM1
               FNM1=FN
               FN=FNP1
   40       CONTINUE
            BESSPH=FN
         ENDIF
      ENDIF
      END
