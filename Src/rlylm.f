C $Id: rlylm.f,v 1.4 1999/01/31 11:43:38 emilio Exp $

      SUBROUTINE RLYLM( LMAX, R, RLY, GRLY )
C *********************************************************************
C FINDS REAL SPHERICAL HARMONICS MULTIPLIED BY R**L: RLY=R**L*YLM,
C AND THEIR GRADIENTS GRLY, AT POINT R:
C    YLM = C * PLM( COS(THETA) ) * SIN(M*PHI)   FOR   M <  0
C    YLM = C * PLM( COS(THETA) ) * COS(M*PHI)   FOR   M >= 0
C WITH (THETA,PHI) THE POLAR ANGLES OF R, C A POSITIVE NORMALIZATION
C CONSTANT AND PLM ASSOCIATED LEGENDRE POLYNOMIALS.
C THE ORDER OF THE Y'S IS THAT IMPLICIT IN THE NESTED LOOPS
C    DO L=0,LMAX
C      DO M=-L,L
C WITH A UNIFIED INDEX ILM=1,...,LMAX**2 INCREASING BY ONE UNIT IN THE
C  INNER LOOP.
C WRITTEN BY J.M.SOLER. AUG/96
C *********** INPUT ***************************************************
C INTEGER LMAX : Maximum angular momentum quantum number required.
C REAL*8  R(3) : Position at which Y and GY are required.
C *********** OUTPUT **************************************************
C REAL*8 RLY(LMAX*LMAX)    : Real spherical harmonics times r**l,
C                             at point R, as explained above.
C REAL*8 GRLY(3,LMAX*LMAX) : Gradient of the RLY's at point R.
C *********** UNITS ***************************************************
C Units of R are arbitrary. Units of RLY and GRLY are related to those
C  of R in the obvious way.
C *********************************************************************

C Next line is non-standard and may be suppressed
      IMPLICIT NONE

C Declare argument types and dimensions
      INTEGER           LMAX
      DOUBLE PRECISION  R(3), RLY(0:*), GRLY(3,0:*)

C Dimension parameters for internal variables
      INTEGER MAXL, MAXLP1
      PARAMETER ( MAXL = 10, MAXLP1 = MAXL+1 )

C Other internal parameters
      DOUBLE PRECISION TINY, ZERO, HALF, ONE, TWO, THREE, SIX
      PARAMETER (TINY=1.D-14, ZERO=0.D0, HALF=0.5D0, ONE=1.D0,
     .           TWO=2.D0, THREE=3.D0, SIX=6.D0)

C Internal variables
      INTEGER
     .  I, ILM, ILM0, L, LMXMX, M, MS
      DOUBLE PRECISION
     .  C(0:MAXLP1*MAXLP1), COSM, COSMM1, COSPHI,
     .  DP(0:MAXLP1,0:MAXLP1), FAC, FOURPI, GY(3),
     .  P(0:MAXLP1,0:MAXLP1),
     .  RL(-1:MAXL), RSIZE, RX, RY, RZ, RXY,
     .  SINM, SINMM1, SINPHI, Y
      SAVE LMXMX, C
      DATA LMXMX /-1/

C Evaluate normalization constants once and for all
      IF (LMAX.GT.LMXMX) THEN
         IF (LMAX.GT.MAXL) STOP 'YLM: MAXL too small'
         FOURPI=TWO**4*ATAN(ONE)
         DO 20 L=0,LMAX
            ILM0=L*L+L
            DO 15 M=0,L
               FAC=(2*L+1)/FOURPI
               DO 10 I=L-M+1,L+M
                  FAC=FAC/I
   10          CONTINUE
               C(ILM0+M)=SQRT(FAC)
C              Next line because Y's are real combinations of M and -M
               IF (M.NE.0) C(ILM0+M)=C(ILM0+M)*SQRT(TWO)
               C(ILM0-M)=C(ILM0+M)
   15       CONTINUE
   20    CONTINUE
         LMXMX=LMAX
      ENDIF

C Initalize to zero
      DO 25 ILM = 0,(LMAX+1)*(LMAX+1)-1
        RLY(ILM) = ZERO
        GRLY(1,ILM) = ZERO
        GRLY(2,ILM) = ZERO
        GRLY(3,ILM) = ZERO
   25 CONTINUE

C Explicit formulas up to L=2
      IF (LMAX.LE.2) THEN

         RLY(0) = C(0)

C        Label 999 is the exit point
         IF (LMAX.EQ.0) GOTO 999

         RLY(1)    = -(C(1)*R(2))
         GRLY(2,1) = -C(1)

         RLY(2)    =  C(2)*R(3)
         GRLY(3,2) =  C(2)

         RLY(3)    = -(C(3)*R(1))
         GRLY(1,3) = -C(3)

         IF (LMAX.EQ.1) GOTO 999

         RLY(4)    =  C(4)*SIX*R(1)*R(2)
         GRLY(1,4) =  C(4)*SIX*R(2)
         GRLY(2,4) =  C(4)*SIX*R(1)

         RLY(5)    = (-C(5))*THREE*R(2)*R(3)
         GRLY(2,5) = (-C(5))*THREE*R(3)
         GRLY(3,5) = (-C(5))*THREE*R(2)

         RLY(6)    =  C(6)*HALF*(TWO*R(3)*R(3)-R(1)*R(1)-R(2)*R(2))
         GRLY(1,6) = (-C(6))*R(1)
         GRLY(2,6) = (-C(6))*R(2)
         GRLY(3,6) =  C(6)*TWO*R(3)

         RLY(7)    = (-C(7))*THREE*R(1)*R(3)
         GRLY(1,7) = (-C(7))*THREE*R(3)
         GRLY(3,7) = (-C(7))*THREE*R(1)

         RLY(8)    =  C(8)*THREE*(R(1)*R(1)-R(2)*R(2))
         GRLY(1,8) =  C(8)*SIX*R(1)
         GRLY(2,8) = (-C(8))*SIX*R(2)

         GOTO 999
      ENDIF

C Special case for R=0
      RSIZE = SQRT( R(1)*R(1)+R(2)*R(2)+R(3)*R(3) )
      IF ( RSIZE .LT. TINY ) THEN
        RLY(0) = C(0)
        GRLY(2,1) = -C(1)
        GRLY(3,2) =  C(2)
        GRLY(1,3) = -C(3)
        GOTO 999
      ENDIF

C Avoid z axis
      RX = R(1) / RSIZE
      RY = R(2) / RSIZE
      RZ = R(3) / RSIZE
      RXY = SQRT(RX*RX+RY*RY)
      IF (RXY.LT.TINY) THEN
        RX = TINY
        RXY = SQRT(RX*RX+RY*RY)
      ENDIF

C General algorithm based on routine PLGNDR of 'Numerical Recipes'
C See also J.M.Soler notes of 23/04/96.
C     Find associated Legendre polynomials and their derivative
      DO 60 M=LMAX,0,-1
         P(M,M+1)=ZERO
         P(M,M)=ONE
         FAC=ONE
         DO 30 I=1,M
            P(M,M)=-(P(M,M)*FAC*RXY)
            FAC=FAC+TWO
   30    CONTINUE
         P(M+1,M)=RZ*(2*M+1)*P(M,M)
         DO 40 L=M+2,LMAX
            P(L,M)=(RZ*(2*L-1)*P(L-1,M)-(L+M-1)*P(L-2,M))/(L-M)
   40    CONTINUE
         DO 50 L=M,LMAX
            DP(L,M)=-((M*P(L,M)*RZ/RXY+P(L,M+1))/RXY)
   50   CONTINUE
   60 CONTINUE
C     Find spherical harmonics and their gradient
      RL(-1) = ZERO
      RL(0)  = ONE
      DO 70 L = 1,LMAX
        RL(L) = RL(L-1)*RSIZE
   70 CONTINUE
      COSPHI=RX/RXY
      SINPHI=RY/RXY
      COSM=ONE
      SINM=ZERO
      DO 90 M=0,LMAX
        DO 80 L=M,LMAX
          DO 75 MS = -1,1,2
            IF (MS.EQ.-1) THEN
              ILM=L*L+L-M
              Y=C(ILM)*P(L,M)*SINM
              GY(1)=(-DP(L,M))*RX *RZ *SINM - P(L,M)*M*COSM*SINPHI/RXY
              GY(2)=(-DP(L,M))*RY *RZ *SINM + P(L,M)*M*COSM*COSPHI/RXY
              GY(3)= DP(L,M)*RXY*RXY*SINM
            ELSE
              ILM=L*L+L+M
              Y=C(ILM)*P(L,M)*COSM
              GY(1)=(-DP(L,M))*RX *RZ *COSM + P(L,M)*M*SINM*SINPHI/RXY
              GY(2)=(-DP(L,M))*RY *RZ *COSM - P(L,M)*M*SINM*COSPHI/RXY
              GY(3)= DP(L,M)*RXY*RXY*COSM
            ENDIF
            GY(1)= GY(1)*C(ILM)/RSIZE
            GY(2)= GY(2)*C(ILM)/RSIZE
            GY(3)= GY(3)*C(ILM)/RSIZE
            RLY(ILM)=RL(L)*Y
            GRLY(1,ILM)= RX*L*RL(L-1)*Y + RL(L)*GY(1)
            GRLY(2,ILM)= RY*L*RL(L-1)*Y + RL(L)*GY(2)
            GRLY(3,ILM)= RZ*L*RL(L-1)*Y + RL(L)*GY(3)
   75     CONTINUE
   80   CONTINUE
        COSMM1=COSM
        SINMM1=SINM
        COSM=COSMM1*COSPHI-SINMM1*SINPHI
        SINM=COSMM1*SINPHI+SINMM1*COSPHI
   90 CONTINUE

  999 CONTINUE
      END
