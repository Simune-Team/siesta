      LOGICAL FUNCTION PROPOR( N, A, B, TOL, AOVERB )
C **********************************************************************
C Checks if two vectors are proportional within a tolerance.
C Written by J.M.Soler. August 1996.
C **********************************************************************

      IMPLICIT NONE

      INTEGER           I, IMAX, N, Node
      DOUBLE PRECISION  A(N), AOVERB, B(N), BMAX, TOL

      BMAX = 0.D0
      IMAX = 0
      DO 10 I = 1,N
        IF ( ABS(B(I)) .GT. BMAX ) THEN
          IMAX = I
          BMAX = ABS(B(I))
        ENDIF
   10 CONTINUE
      IF (IMAX .EQ. 0) call die("propor: ERROR:  IMAX = 0")

      PROPOR = .TRUE.
      IF (BMAX .EQ. 0.D0) THEN
        AOVERB = 0.D0
        DO 20 I = 1,N
          IF ( ABS(A(I)) .GT. TOL ) THEN
            PROPOR = .FALSE.
            RETURN
          ENDIF
   20   CONTINUE
      ELSE
        AOVERB = A(IMAX) / B(IMAX)
        DO 30 I = 1,N
          IF ( ABS(A(I)-B(I)*AOVERB) .GT. TOL ) THEN
            PROPOR = .FALSE.
            RETURN
          ENDIF
   30   CONTINUE
      ENDIF
      END



