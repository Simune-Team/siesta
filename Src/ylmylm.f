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

      INTEGER MAXLM

      INTEGER           I, L, LOFILM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::
     .  Y
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     .  DYDR

      EXTERNAL          LOFILM, MEMORY

      L = MAX( LOFILM(ILM1), LOFILM(ILM2) )
      MAXLM = (L+1)*(L+1)

      allocate(Y(MAXLM))
      call memory('A','D',MAXLM,'ylmylm')
      allocate(DYDR(3,MAXLM))
      call memory('A','D',3*MAXLM,'ylmylm')

      CALL RLYLM( L, R, Y, DYDR )

      YY = Y(ILM1) * Y(ILM2)
      DO 10 I = 1,3
        DYYDR(I) = DYDR(I,ILM1) * Y(ILM2) + Y(ILM1) * DYDR(I,ILM2)
   10 CONTINUE

      call memory('D','D',size(Y),'ylmylm')
      deallocate(Y)
      call memory('D','D',size(DYDR),'ylmylm')
      deallocate(DYDR)

      END

