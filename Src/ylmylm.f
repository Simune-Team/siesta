      SUBROUTINE YLMYLM( ILM1, ILM2, R, YY, DYYDR )

C *********************************************************************
C Returns the product of two real spherical harmonics (SH),
C  times r**l each, and its gradient.
C *********** INPUT ***************************************************
C INTEGER ILM1, ILM2 : Combined SH indexes. ILM = L*L+L+M+1
C REAL*8  R(3)       : Vector towards the (theta,phi) direction.
C *********** OUTPUT **************************************************
C REAL*8  YY       : Product of the two SH.
C REAL*8  DYYDR(3) : Derivative (gradient) of YY with respect to R
C *********************************************************************
C Written by J.M.Soler. Feb' 96.
C *********************************************************************

      IMPLICIT          NONE
      INTEGER           ILM1, ILM2
      DOUBLE PRECISION  DYYDR(3), R(3), YY

      INTEGER MAXL
      PARAMETER (MAXL = 6)

      INTEGER           I, L, LOFILM
      DOUBLE PRECISION  DYDR(3,MAXL*MAXL), Y(MAXL*MAXL)
      EXTERNAL          CHKDIM, LOFILM

      L = MAX( LOFILM(ILM1), LOFILM(ILM2) )
      CALL CHKDIM( 'YLMYLM', 'MAXL', MAXL, L, 1 )
      CALL RLYLM( L, R, Y, DYDR )

      YY = Y(ILM1) * Y(ILM2)
      DO 10 I = 1,3
        DYYDR(I) = DYDR(I,ILM1) * Y(ILM2) + Y(ILM1) * DYDR(I,ILM2)
   10 CONTINUE

      END

