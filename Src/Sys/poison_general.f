C $Id: poison_general.f,v 1.1 1999/02/26 10:28:10 wdpgaara Exp $

      SUBROUTINE POISON( CELL, N1, N2, N3, RHO, U, V, STRESS, NCG, CG )

C *********************************************************************
C SOLVES POISSON'S EQUATION.
C ENERGY AND POTENTIAL RETURNED IN RYDBERG UNITS.
C WRITTEN BY J.M.SOLER, JUNE 1995.
C **************** INPUT **********************************************
C REAL*8  CELL(3,3)     : UNIT CELL VECTORS
C INTEGER N1,N2,N3      : NUMBER OF MESH DIVISIONS IN EACH CELL VECTOR
C REAL*4  RHO(N1,N2,N3) : DENSITIY AT MESH POINTS
C **************** OUTPUT *********************************************
C REAL*8  U             : ELECTROSTATIC ENERGY (IN RY)
C REAL*4  V(N1,N2,N3)   : ELECTROSTATIC POTENTIAL (IN RY)
C                         V AND RHO MAY BE THE SAME PHYSICAL ARRAY
C REAL*8  STRESS(3,3) : Electrostatic-energy contribution to stress
C                       tensor (in Ry/Bohr**3) assuming constant density
C                       (not charge), i.e. r->r' => rho'(r') = rho(r)
C                       For plane-wave and grid (finite difference)
C                       basis sets, density rescaling gives an extra
C                       term (not included) equal to -2*U/cell_volume
C                       for the diagonal elements of stress. For other
C                       basis sets, the extra term is, in general:
C                       IntegralOf( V * d_rho/d_strain ) / cell_volume
C **************** INPUT AND OUTPUT ***********************************
C INTEGER NCG : Dimension of array CG. It must be 2*N1*N2*N3
C               or larger. If it is smaller, no calculations are
C               performed, and it is set on output to this value.
C               If it is large enough on input, it is not changed.
C **************** AUXILIARY ******************************************
C REAL*8  CG(NCG) : SPACE THAT CAN BE USED FREELY OUTSIDE.
C *********************************************************************

      IMPLICIT          NONE
      INTEGER           N1, N2, N3, NCG
      REAL              RHO(N1*N2*N3), V(N1*N2*N3)
      DOUBLE PRECISION  CELL(3,3), STRESS(3,3), U
      DOUBLE PRECISION  CG(2,N1*N2*N3)

      INTEGER           I, I1, I2, I3, IX, J, J1, J2, J3, JX,
     .                  MESH(3), NP
      DOUBLE PRECISION  C, B(3,3), DU, G(3), G2, G2MAX, 
     .                  K0(3), PI, TINY, VG, VOLUME, VOLCEL
      EXTERNAL          CHKGMX, FFT, RECLAT, VOLCEL
      SAVE              K0, TINY

      DATA K0   /3*0.D0/
      DATA TINY /1.D-15/

C     START TIME COUNTER
      CALL TIMER ('POISON',1)

C     CHECK THAT CG IS LARGE ENOUGH
      IF (NCG .LT. 2*N1*N2*N3) THEN
        NCG = 2*N1*N2*N3
C       GO TO EXIT POINT
        GOTO 999
      ENDIF

C     FIND UNIT CELL VOLUME
      VOLUME = VOLCEL( CELL )

C     FIND RECIPROCAL LATTICE VECTORS
      CALL RECLAT (CELL, B, 1 )

C     FIND MAXIMUM PLANEWAVE CUTOFF
      NP = N1 * N2 * N3
      MESH(1) = N1
      MESH(2) = N2
      MESH(3) = N3
      G2MAX = 1.D30
      CALL CHKGMX( K0, B, MESH, G2MAX )

C     COPY DENSITY TO COMPLEX ARRAY
      DO 5 I = 1, NP
        CG(1,I) = RHO(I)
        CG(2,I) = 0.D0
    5 CONTINUE

C     FIND FOURIER TRANSFORM OF DENSITY
      CALL FFT( CG, MESH, -1 )
 
C     INITIALIZE STRESS CONTRIBUTION
      DO 15 IX = 1,3
        DO 10 JX = 1,3
          STRESS(JX,IX) = 0.D0
   10   CONTINUE
   15 CONTINUE

C     MULTIPLY BY 8*PI/G2 TO GET THE POTENTIAL
      PI = 4.D0 * ATAN(1.D0)
      U = 0.D0
      DO 40 I3 = -((N3-1)/2), N3/2
      DO 40 I2 = -((N2-1)/2), N2/2
      DO 40 I1 = -((N1-1)/2), N1/2
        G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
        G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
        G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3
        G2 = G(1)**2 + G(2)**2 + G(3)**2
        J1 = MOD( I1 + N1, N1 )
        J2 = MOD( I2 + N2, N2 )
        J3 = MOD( I3 + N3, N3 )
        J = 1 + J1 + N1 * J2 + N1 * N2 * J3
        IF (G2.LT.G2MAX .AND. G2.GT.TINY) THEN
          VG = 8.D0 * PI / G2
          DU = VG * ( CG(1,J)**2 + CG(2,J)**2 )
          U = U + DU
          C = 2.D0 * DU / G2
          DO 30 IX = 1,3
            DO 20 JX = 1,3
              STRESS(JX,IX) = STRESS(JX,IX) + C * G(IX) * G(JX)
   20       CONTINUE
   30     CONTINUE
          CG(1,J) = VG * CG(1,J)
          CG(2,J) = VG * CG(2,J)
        ELSE
          CG(1,J) = 0.D0
          CG(2,J) = 0.D0
        ENDIF
   40 CONTINUE
      U = 0.5D0 * U * VOLUME / DBLE(N1*N2*N3)**2
      C = 0.5D0 / DBLE(N1*N2*N3)**2
      DO 60 IX = 1,3
        DO 50 JX = 1,3
          STRESS(JX,IX) = C * STRESS(JX,IX)
   50   CONTINUE
        STRESS(IX,IX) = STRESS(IX,IX) + U / VOLUME
   60 CONTINUE
 
C     GO BACK TO REAL SPACE
      CALL FFT( CG, MESH, +1 )

C     COPY POTENTIAL TO ARRAY V
      DO 70 I = 1, NP
        V(I) = CG(1,I)
   70 CONTINUE
 
C     STOP TIME COUNTER
  999 CONTINUE
      CALL TIMER ('POISON',2)
      END



      SUBROUTINE FFT (F,MESH,ISN)
C        THIS IS AN INTERFACE TO ROUTINE CFT 
C        Notice that ISN=1 means INVERSE trnasform and that
C        the factor 1/N is in the inverse transform.
C        Written by J.M.Soler. Last revision April 1997.

C        FOR SINGLE PRECISION, SIMPLY COMMENT OUT NEXT LINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(*),MESH(3)
      PARAMETER (ONE=1.D0)
      CALL TIMER ('FFT   ',1)
      N1=MESH(1)
      N2=MESH(2)
      N3=MESH(3)
      N=N1*N2*N3
      CALL CFT (F,F(2), N, N1, N1,       2*ISN)
      CALL CFT (F,F(2), N, N2, N1*N2,    2*ISN)
      CALL CFT (F,F(2), N, N3, N1*N2*N3, 2*ISN)
*     IF (ISN.LT.0) THEN
      IF (ISN.GT.0) THEN
        C=ONE/N
        DO 10 I=1,2*N
          F(I)=F(I)*C
   10   CONTINUE
      ENDIF
      CALL TIMER ('FFT   ',2)
      END



      SUBROUTINE NFFT( N )

C CHANGES N INTO THE NEXT INTEGER ALLOWED BY THE FFT-ROUTINE CFT
C WRITTEN BY J.M.SOLER. MAY/95.

      PARAMETER (NP = 3, NMAX = 1000000)
      INTEGER IPRIME(NP)
      DATA IPRIME / 2, 3, 5 /

      NMIN = N
      DO 30 N = NMIN, NMAX
        NREM = N
        DO 20 IP = 1,NP
   10     CONTINUE 
          IF ( MOD( NREM, IPRIME(IP) ) .EQ. 0 ) THEN
            NREM = NREM / IPRIME(IP)
            GOTO 10
          ENDIF
   20   CONTINUE
        IF (NREM .EQ. 1) RETURN
   30 CONTINUE
      WRITE(6,*) 'NFFT: NO SUITABLE INTEGER FOUND FOR N =', NMIN
      STOP
      END
