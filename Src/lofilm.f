      INTEGER FUNCTION LOFILM( ILM )

C *********************************************************************
C Finds the total angular momentum quantum number L which corresponds
C to the combined index ILM == (L,M) implicit in the nested loop
C    ILM = 0
C    DO L=0,LMAX
C      DO M=-L,L
C        ILM = ILM + 1
C with ILM = 1,2,...,L**2
C Written by J.M.Soler. April 1996.
C *********************************************************************

      IMPLICIT NONE

      INTEGER ILM, JLM, MAXL
      PARAMETER ( MAXL = 100 )

      IF ( ILM .LE. 0 ) STOP 'LOFILM: ILM not allowed'
      JLM = 0
      DO 10 LOFILM = 0,MAXL
        JLM = JLM + 2*LOFILM + 1
        IF ( JLM .GE. ILM ) RETURN
   10 CONTINUE
      STOP 'LOFILM: ILM too large'
      END

