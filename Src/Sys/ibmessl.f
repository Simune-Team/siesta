      include '../rdiag.f'


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
C REAL*8  U           : ELECTROSTATIC ENERGY (IN RY)
C REAL*4  V(N1,N2,N3) : ELECTROSTATIC POTENTIAL (IN RY)
C                       V AND RHO MAY BE THE SAME PHYSICAL ARRAY
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
C INTEGER NCG : Dimension of array CG. It must be (N1/2+1)*N2*N3
C               or larger. If it is smaller, no calculations are
C               performed, and it is set on output to this value.
C               If it is large enough on input, it is not changed.
C **************** AUXILIARY ******************************************
C REAL*8  CG(NCG) : SPACE THAT CAN BE USED FREELY OUTSIDE.
C                   AUX MAY BE THE SAME PHYSICAL ARRAY AS RHO
C                   AND/OR V BUT NOTICE THAT IT IS LARGER
C *********************************************************************

      IMPLICIT          NONE
      INTEGER           N1, N2, N3, NCG
      REAL              RHO(N1*N2*N3), V(N1*N2*N3)
      DOUBLE PRECISION  CELL(3,3), STRESS(3,3), U
      COMPLEX           CG( 0:N1/2, 0:N2-1, 0:N3-1 )

      INTEGER           I1, I2, I3, IX, J1, J2, J3, JX,
     .                  MESH(3)
      DOUBLE PRECISION  C, B(3,3), DU, G(3), G2, G2MAX, 
     .                  K0(3), PI, TINY, VG, VOLUME, VOLCEL
      EXTERNAL          CHKGMX, RFFT, RECLAT, VOLCEL
      SAVE              K0, TINY

      DATA K0   /3*0.D0/
      DATA TINY /1.D-15/

C     START TIME COUNTER
      CALL TIMER ('POISON',1)

C     CHECK THAT CG IS LARGE ENOUGH
      IF (NCG .LT. (N1/2+1)*N2*N3) THEN
        NCG = (N1/2+1)*N2*N3
C       GO TO EXIT POINT
        GOTO 999
      ENDIF

C     FIND UNIT CELL VOLUME
      VOLUME = VOLCEL( CELL )

C     FIND RECIPROCAL LATTICE VECTORS
      CALL RECLAT (CELL, B, 1 )

C     FIND MAXIMUN PLANEWAVE CUTOFF
      MESH(1) = N1
      MESH(2) = N2
      MESH(3) = N3
      G2MAX = 1.D30
      CALL CHKGMX( K0, B, MESH, G2MAX )

C     FIND FOURIER TRANSFORM OF DENSITY
      CALL RFFT( RHO, N1, N2, N3, CG, +1 )
 
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
      DO 40 I1 = 0, N1/2
        G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
        G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
        G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3
        G2 = G(1)**2 + G(2)**2 + G(3)**2
        J1 = MOD( I1 + N1, N1 )
        J2 = MOD( I2 + N2, N2 )
        J3 = MOD( I3 + N3, N3 )
        IF (G2.LT.G2MAX .AND. G2.GT.TINY) THEN
          VG = 8.D0 * PI / G2
          DU = VG * CABS( CG(J1,J2,J3) )**2
          IF (J1.NE.0 .AND. J1.NE.N1/2) DU = 2.D0 * DU
          U = U + DU
          C = 2.D0 * DU / G2
          DO 30 IX = 1,3
            DO 20 JX = 1,3
              STRESS(JX,IX) = STRESS(JX,IX) + C * G(IX) * G(JX)
   20       CONTINUE
   30     CONTINUE
          CG(J1,J2,J3) = VG * CG(J1,J2,J3)
        ELSE
          CG(J1,J2,J3) = (0.D0, 0.D0)
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
      CALL RFFT( V, N1, N2, N3, CG, -1 )
 
C     STOP TIME COUNTER
  999 CONTINUE
      CALL TIMER ('POISON',2)
      END



      SUBROUTINE RFFT( R, N1, N2, N3, C, ISIG )

C **********************************************************************
C PERFORMS A SINGLE-PRECISION 3D FFT OF A REAL FUNCTION
C THIS IS AN INTERFACE TO ROUTINES SRCFT3 AND SCRFT3 OF IBM-ESSL LIBRARY
C Notice that the factor 1/N is in the inverse transform.
C WRITTEN BY J.M.SOLER. Last revision April 1997.
C ***** INPUT **********************************************************
C INTEGER N1,N2,N3 : NUMBER OF POINTS IN EACH DIRECTION
C INTEGER ISIG     : FOURIER-TRANSFORM DIRECTION
C                    ISIG = +1  =>  FROM REAL TO FOURIER SPACE
C                    ISIG = -1  =>  FROM FOURIER TO REAL SPACE
C ***** INPUT OR OUTPUT (DEPENDING OF ISIG) ****************************
C REAL*4    R(0:N1-1,0:N2-1,0:N3-1) : REAL FUNCTION VALUES
C COMPLEX*8 C(0:N1/2,0:N2-1,0:N3-1) : COMPLEX FOURIER COEFFICIENTS
C                               R AND C MAY BE THE SAME PHYSICAL ARRAY,
C                               BUT NOTICE THAT C IS LARGER
C **********************************************************************

      IMPLICIT   NONE
      INTEGER    N1, N2, N3, ISIG
      REAL       R(0:N1-1,0:N2-1,0:N3-1)
      COMPLEX    C( 0:N1/2, 0:N2-1, 0:N3-1 )

      INTEGER    NAUX, INC2R, INC3R, INC2C, INC3C
      REAL       ONE
      PARAMETER (NAUX=65000, ONE=1.D0)
      REAL       AUX(NAUX,7), SCALE

      CALL TIMER ('RFFT',1)

      INC2R = N1 + 2
      INC3R = INC2R * N2
      INC2C = INC2R / 2
      INC3C = INC3R / 2

      IF ( ISIG .GT. 0 ) THEN
C        REAL TO FOURIER SPACE
         SCALE = ONE
         CALL COPYM( N1, N2*N3, R, N1+2, N2*N3, C )
         CALL SRCFT3( C, INC2R, INC3R, C, INC2C, INC3C, N1, N2, N3,
     .                ISIG, SCALE, AUX, NAUX )
      ELSE
C        FOURIER TO REAL SPACE
         SCALE = ONE / (N1*N2*N3)
         CALL SCRFT3( C, INC2C, INC3C, C, INC2R, INC3R, N1, N2, N3,
     .                ISIG, SCALE, AUX, NAUX )
         CALL COPYM( N1+2, N2*N3, C, N1, N2*N3, R )
      ENDIF

      CALL TIMER ('RFFT',2)
      END



      SUBROUTINE COPYM( N1A, N2A, A, N1B, N2B, B )

C COPIES MATRIX A(N1A,N2A) INTO MATRIX B(N1B,N2B).
C A AND B MAY BE THE SAME PHYSICAL ARRAY.
C WRITTEN BY J.M.SOLER. JUNE/95.

      IMPLICIT   NONE
      INTEGER    N1A, N2A, N1B, N2B
      REAL       A(N1A,N2A), B(N1B,N2B)

      INTEGER    I1, I2
      REAL       ZERO
      PARAMETER  (ZERO = 0.D0)

      IF ( N1B .GT. N1A ) THEN
        DO 20 I2 = N2B, N2A+1, -1
          DO 10 I1 = N1B, 1, -1
            B(I1,I2) = ZERO
   10     CONTINUE
   20   CONTINUE
        DO 50 I2 = N2A, 1, -1
          DO 30 I1 = N1B, N1A+1, -1
            B(I1,I2) = ZERO
   30     CONTINUE
          DO 40 I1 = N1A, 1, -1
            B(I1,I2) = A(I1,I2)
   40     CONTINUE
   50   CONTINUE
      ELSE
        DO 70 I2 = 1, N2A
          DO 60 I1 = 1, N1B
            B(I1,I2) = A(I1,I2)
   60     CONTINUE
   70   CONTINUE
        DO 90 I2 = N2A+1, N2B
          DO 80 I1 = 1, N1B
            B(I1,I2) = ZERO
   80     CONTINUE
   90   CONTINUE
      ENDIF

      END



      SUBROUTINE NFFT( N )

C CHANGES N INTO THE NEXT INTEGER ALLOWED BY THE FFT ROUTINE
C THIS IS A PROVISIONAL VERSION FOR VALUES OF N UP TO 2048.
C REF: IBM-ESSL LIBRARY MANUAL.
C WRITTEN BY J.M.SOLER. JUNE/95.

      PARAMETER (NNDATA = 63)
      INTEGER N, NDATA(NNDATA)
      DATA NDATA / 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 30, 32, 36, 40,
     .  42, 48, 56, 60, 64, 70, 72, 80, 84, 90, 96, 112, 120, 126, 128,
     .  140, 144, 160, 168, 180, 192, 210, 224, 240, 252, 256, 280, 288,
     .  320, 336, 360, 384, 420, 448, 480, 504, 512, 560, 576, 630,
     .  640, 672, 720, 768, 840, 896, 960, 1008, 1024 /

      DO 10 IND = 1, NNDATA
        IF (2*NDATA(IND) .GE. N) THEN
          N = 2*NDATA(IND)
          RETURN
        ENDIF
   10 CONTINUE

      WRITE(6,*) 'NFFT: TOO LARGE N =', N
      STOP
      END

      SUBROUTINE CPUTIM (TIME)

C  RETURNS CPU TIME IN SECONDS SINCE PROGRAM START
C  WRITEN BY J.SOLER (JSOLER AT EMDUAM11)

      DOUBLE PRECISION TIME
*     REAL TIMES(2)

C  NEXT LINES FOR IBM SYSTEMS-370 (ASSEMBLE ROUTINE TASKTM REQUIRED)
*     CALL TASKTM (ITIME)
*     TIME = 1.D-4 * ITIME

C  NEXT LINES FOR IBM-3090
*     CALL CPUTIME (TIME,RCODE)
*     TIME = 1.D-6 * TIME

C  NEXT LINE FOR CRAY
*     TIME = SECOND()

C  NEXT LINE FOR SUN OR DEC WORKSTATIONS
*     TIME = ETIME(TIMES)

C  NEXT LINE FOR IBM RS/6000 WORKSTATIONS
      TIME = MCLOCK()*0.01D0

C  NEXT LINE FOR ANYTHING ELSE
*     TIME = 0.D0

      END
