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

C     FIND MAXIMUN PLANEWAVE CUTOFF
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
C        THIS IS AN INTERFACE TO ROUTINE CFT BELOW
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

      SUBROUTINE CFT(A,B,NTOT,N,NSPAN,ISN)
C
C     MULTIVARIATE COMPLEX FOURIER TRANSFORM, COMPUTED IN PLACE
C     USING MIXED-RADIX FAST FOURIER TRANSFORM ALGORITHM.
C     BY R. C. SINGLETON, STANFORD RESEARCH INSTITUTE, OCT. 1968
C     ARRAYS A AND B ORIGINALLY HOLD THE REAL AND IMAGINARY
C     COMPONENTS OF THE DATA, AND RETURN THE REAL AND
C     IMAGINARY COMPONENTS OF THE RESULTING FOURIER COEFFICIENTS.
C     MULTIVARIATE DATA IS INDEXED ACCORDING TO THE FORTRAN
C     ARRAY ELEMENT SUCCESSOR FUNCTION, WITHOUT LIMIT
C     ON THE NUMBER OF IMPLIED MULTIPLE SUBSCRIPTS.
C     THE SUBROUTINE IS CALLED ONCE FOR EACH VARIATE.
C     THE CALLS FOR A MULTIVARIATE TRANSFORM MAY BE IN ANY ORDER.
C     NTOT IS THE TOTAL NUMBER OF COMPLEX DATA VALUES.
C     N IS THE DIMENSION OF THE CURRENT VARIABLE.
C     NSPAN/N IS THE SPACING OF CONSUCUTIVE DATA VALUES
C     WHILE INDEXING THE CURRENT VARIABLE.
C     THE SIGN OF ISN DETERMINES THE SIGN OF THE COMPLEX
C     EXPONENTIAL, AND THE MAGNITUDE OF ISN IS NORMALLY ONE.
C
C     FOR A SINGLE-VARIATE TRANSFORM,
C     NTOT = N = NSPAN = (NUMBER OF COMPLEX DATA VALUES), F.G.
C     CALL CFT(A,B,N,N,N,1)
C
C     A TRI-VARIATE TRANSFORM WITH A(N1,N2,N3), B(N1,N2,N3)
C     IS COMPUTED BY
C     CALL CFT(A,B,N1*N2*N3,N1,N1,1)
C     CALL CFT(A,B,N1*N2*N3,N2,N1*N2,1)
C     CALL CFT(A,B,N1*N2*N3,N3,N1*N2*N3,1)
C
C     THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
C     ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
C     GIVE THE CORRECT INDEXING INCREMENT AND THE SECOND PARAMETER
C     USED TO PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF
C     IMAGINARY VALUES, E.G.
C
C        REAL S(2)
C        EQUIVALENCE (A,S)
C        ....
C        ....
C        CALL CFT(A,S(2),NTOT,N,NSPAN,2)
C
C     ARRAYS AT(MAXF), CK(MAXF), BT(MAXF), SK(MAXF), AND NP(MAXP)
C     ARE USED FOR TEMPORARY STORAGE. IF THE AVAILABLE STORAGE
C     IS INSUFFICIENT, THE PROGRAM IS TERMINATED BY A STOP.
C     MAXF MUST BE .GE. THE MAXIMUM PRIME FACTOR OF N.
C     MAXP MUST BE .GT. THE NUMBER OF PRIME FACTORS OF N.
C     IN ADDITION, IF THE SQUARE-FREE PORTION K CF N HAS TWO OR
C     MORE PRIME FACTORS, THEN MAXP MUST BE .GE. K-1.
C     ARRAY STORAGE IN NFAC FOR A MAXIMUM OF 11 FACTORS OF N.
C     IF N HAS MORE THAN ONE SQUARE-FREE FACTOR, THE PRODUCT OF THE
C     SQUARE-FREE FACTORS MUST BE .LE. 210

*        ADAPTED TO DOUBLE PRECISION BY J.SOLER 30/10/89.
*        FOR SINGLE PRECISION SIMPLY COMMENT OUT NEXT LINE.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*        PARAMETERS INTRODUCED BY JSOLER
C        ARRAY STORAGE FOR MAXIMUM PRIME FACTOR OF 23
      PARAMETER (MXFDIM=23,MXPDIM=209)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0)
      DIMENSION A(*),B(*)
      DIMENSION NFAC(11),NP(MXPDIM)
      DIMENSION AT(MXFDIM),CK(MXFDIM),BT(MXFDIM),SK(MXFDIM)
*        NEXT LINE SUPRESSED BY J.SOLER AND II CHANGED
*        TO I EVERYWHERE (ORIGINAL LINES COMMENTED OUT)
*        EQUIVALENCE (I,II)
      MAXF=MXFDIM
      MAXP=MXPDIM
      IF(N .LT. 2) RETURN
      INC=ISN
C     THE FOLLOWING CONSTANTS ARE RAD = 2.*PI , S72 = SIN(0.4*PI) ,
C     C72 = COS(0.4*PI) AND S120 = SQRT(0.75)
      RAD = 6.2831853071796D0
      S72 = 0.95105651629515D0
      C72 = 0.30901699437495D0
      S120 = 0.86602540378444D0
      IF(ISN .GE. 0) GO TO 10
      S72=-S72
      S120=-S120
      RAD=-RAD
      INC=-INC
   10 NT=INC*NTOT
      KS=INC*NSPAN
      KSPAN=KS
      NN=NT-INC
      JC=KS/N
      RADF=RAD*JC*HALF
      I=0
      JF=0
C     DETERMINE THE FACTORS OF N
      M=0
      K=N
      GO TO 20
   15 M=M+1
      NFAC(M)=4
      K=K/16
   20 IF(K-(K/16)*16 .EQ. 0) GO TO 15
      J=3
      JJ=9
      GO TO 30
   25 M=M+1
      NFAC(M)=J
      K=K/JJ
   30 IF(MOD(K,JJ) .EQ. 0) GO TO 25
      J=J+2
      JJ=J**2
      IF(JJ .LE. K) GO TO 30
      IF(K .GT. 4) GO TO 40
      KT=M
      NFAC(M+1)=K
      IF(K .NE. 1) M=M+1
      GO TO 80
   40 IF(K-(K/4)*4 .NE. 0) GO TO 50
      M=M+1
      NFAC(M)=2
      K=K/4
   50 KT=M
      J=2
   60 IF(MOD(K,J) .NE. 0) GO TO 70
      M=M+1
      NFAC(M)=J
      K=K/J
   70 J=((J+1)/2)*2+1
      IF(J .LE. K) GO TO 60
   80 IF(KT .EQ. 0) GO TO 100
      J=KT
   90 M=M+1
      NFAC(M)=NFAC(J)
      J=J-1
      IF(J .NE. 0) GO TO 90
C     COMPUTE FOURIER TRANSFORM
  100 SD=RADF/KSPAN
      CD=TWO*SIN(SD)**2
      SD=SIN(SD+SD)
      KK=1
      I=I+1
      IF(NFAC(I) .NE. 2) GO TO 400
C     TRANSFORM FOR FACTOR OF 2 (INCLUDING ROTATION FACTOR)
      KSPAN=KSPAN/2
      K1=KSPAN+2
  210 K2=KK+KSPAN
      AK=A(K2)
      BK=B(K2)
      A(K2)=A(KK)-AK
      B(K2)=B(KK)-BK
      A(KK)=A(KK)+AK
      B(KK)=B(KK)+BK
      KK=K2+KSPAN
      IF(KK .LE. NN) GO TO 210
      KK=KK-NN
      IF(KK .LE. JC) GO TO 210
      IF(KK .GT. KSPAN) GO TO 800
  220 C1=ONE-CD
      S1=SD
  230 K2=KK+KSPAN
      AK=A(KK)-A(K2)
      BK=B(KK)-B(K2)
      A(KK)=A(KK)+A(K2)
      B(KK)=B(KK)+B(K2)
      A(K2)=C1*AK-S1*BK
      B(K2)=S1*AK+C1*BK
      KK=K2+KSPAN
      IF(KK .LT. NT) GO TO 230
      K2=KK-NT
      C1=-C1
      KK=K1-K2
      IF(KK .GT. K2) GO TO 230
      AK=C1-(CD*C1+SD*S1)
      S1=(SD*C1-CD*S1)+S1
C     THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C     ERROR. IF ROUNDED ARITHMETIC IS USED, THEY MAY BE DELETED.
C     C1=HALF/(AK**2+S1**2)+HALF
C     S1=C1*S1
C     C1=C1*AK
C     NEXT STATEMENT SHOULD BE DELETED IF NON-ROUNDED ARITHMETIC IS USED
      C1=AK
      KK=KK+JC
      IF(KK .LT. K2) GO TO 230
      K1=K1+INC+INC
      KK=(K1-KSPAN)/2+JC
      IF(KK .LE. JC+JC) GO TO 220
      GO TO 100
C     TRANSFORM FOR FACTOR OF 3 (OPTIONAL CODE)
  320 K1=KK+KSPAN
      K2=K1+KSPAN
      AK=A(KK)
      BK=B(KK)
      AJ=A(K1)+A(K2)
      BJ=B(K1)+B(K2)
      A(KK)=AK+AJ
      B(KK)=BK+BJ
      AK=(-HALF)*AJ+AK
      BK=(-HALF)*BJ+BK
      AJ=(A(K1)-A(K2))*S120
      BJ=(B(K1)-B(K2))*S120
      A(K1)=AK-BJ
      B(K1)=BK+AJ
      A(K2)=AK+BJ
      B(K2)=BK-AJ
      KK=K2+KSPAN
      IF(KK .LT. NN) GO TO 320
      KK=KK-NN
      IF(KK .LE. KSPAN) GO TO 320
      GO TO 700
C     TRANSFORM FOR FACTOR OF 4
  400 IF(NFAC(I) .NE. 4) GO TO 600
      KSPNN=KSPAN
      KSPAN=KSPAN/4
  410 C1=ONE
      S1=0
  420 K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
      AKP=A(KK)+A(K2)
      AKM=A(KK)-A(K2)
      AJP=A(K1)+A(K3)
      AJM=A(K1)-A(K3)
      A(KK)=AKP+AJP
      AJP=AKP-AJP
      BKP=B(KK)+B(K2)
      BKM=B(KK)-B(K2)
      BJP=B(K1)+B(K3)
      BJM=B(K1)-B(K3)
      B(KK)=BKP+BJP
      BJP=BKP-BJP
      IF(ISN .LT. 0) GO TO 450
      AKP=AKM-BJM
      AKM=AKM+BJM
      BKP=BKM+AJM
      BKM=BKM-AJM
      IF(S1 .EQ. ZERO) GO TO 460
  430 A(K1)=AKP*C1-BKP*S1
      B(K1)=AKP*S1+BKP*C1
      A(K2)=AJP*C2-BJP*S2
      B(K2)=AJP*S2+BJP*C2
      A(K3)=AKM*C3-BKM*S3
      B(K3)=AKM*S3+BKM*C3
      KK=K3+KSPAN
      IF(KK .LE. NT) GO TO 420
  440 C2=C1-(CD*C1+SD*S1)
      S1=(SD*C1-CD*S1)+S1
C     THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C     ERROR. IF ROUNDED ARITHMETIC IS USED, THEY MAY BE DELETED.
C     C1=HALF/(C2**2+S1**2)+HALF
C     S1=C1*S1
C     C1=C1*C2
C     NEXT STATEMENT SHOULD BE DELETED IF NON-ROUNDED ARITHMETIC IS USED
      C1=C2
      C2=C1**2-S1**2
      S2=TWO*C1*S1
      C3=C2*C1-S2*S1
      S3=C2*S1+S2*C1
      KK=KK-NT+JC
      IF(KK .LE. KSPAN) GO TO 420
      KK=KK-KSPAN+INC
      IF(KK .LE. JC) GO TO 410
      IF(KSPAN .EQ. JC) GO TO 800
      GO TO 100
  450 AKP=AKM+BJM
      AKM=AKM-BJM
      BKP=BKM-AJM
      BKM=BKM+AJM
      IF(S1 .NE. ZERO) GO TO 430
  460 A(K1)=AKP
      B(K1)=BKP
      A(K2)=AJP
      B(K2)=BJP
      A(K3)=AKM
      B(K3)=BKM
      KK=K3+KSPAN
      IF(KK .LE. NT) GO TO 420
      GO TO 440
C     TRANSFORM FOR FACTOR OF 5 (OPTIONAL CODE)
  510 C2=C72**2-S72**2
      S2=TWO*C72*S72
  520 K1=KK+KSPAN
      K2=K1+KSPAN
      K3=K2+KSPAN
      K4=K3+KSPAN
      AKP=A(K1)+A(K4)
      AKM=A(K1)-A(K4)
      BKP=B(K1)+B(K4)
      BKM=B(K1)-B(K4)
      AJP=A(K2)+A(K3)
      AJM=A(K2)-A(K3)
      BJP=B(K2)+B(K3)
      BJM=B(K2)-B(K3)
      AA=A(KK)
      BB=B(KK)
      A(KK)=AA+AKP+AJP
      B(KK)=BB+BKP+BJP
      AK=AKP*C72+AJP*C2+AA
      BK=BKP*C72+BJP*C2+BB
      AJ=AKM*S72+AJM*S2
      BJ=BKM*S72+BJM*S2
      A(K1)=AK-BJ
      A(K4)=AK+BJ
      B(K1)=BK+AJ
      B(K4)=BK-AJ
      AK=AKP*C2+AJP*C72+AA
      BK=BKP*C2+BJP*C72+BB
      AJ=AKM*S2-AJM*S72
      BJ=BKM*S2-BJM*S72
      A(K2)=AK-BJ
      A(K3)=AK+BJ
      B(K2)=BK+AJ
      B(K3)=BK-AJ
      KK=K4+KSPAN
      IF(KK .LT. NN) GO TO 520
      KK=KK-NN
      IF(KK .LE. KSPAN) GO TO 520
      GO TO 700
C     TRANSFORM FOR ODD FACTORS
  600 K=NFAC(I)
      KSPNN=KSPAN
      KSPAN=KSPAN/K
      IF(K .EQ. 3) GO TO 320
      IF(K .EQ. 5) GO TO 510
      IF(K .EQ. JF) GO TO 640
      JF=K
      S1=RAD/K
      C1=COS(S1)
      S1=SIN(S1)
      IF(JF .GT. MAXF) GO TO 998
      CK(JF)=ONE
      SK(JF)=ZERO
      J=1
  630 CK(J)=CK(K)*C1+SK(K)*S1
      SK(J)=CK(K)*S1-SK(K)*C1
      K=K-1
      CK(K)=CK(J)
      SK(K)=-SK(J)
      J=J+1
      IF(J .LT. K) GO TO 630
  640 K1=KK
      K2=KK+KSPNN
      AA=A(KK)
      BB=B(KK)
      AK=AA
      BK=BB
      J=1
      K1=K1+KSPAN
  650 K2=K2-KSPAN
      J=J+1
      AT(J)=A(K1)+A(K2)
      AK=AT(J)+AK
      BT(J)=B(K1)+B(K2)
      BK=BT(J)+BK
      J=J+1
      AT(J)=A(K1)-A(K2)
      BT(J)=B(K1)-B(K2)
      K1=K1+KSPAN
      IF(K1 .LT. K2) GO TO 650
      A(KK)=AK
      B(KK)=BK
      K1=KK
      K2=KK+KSPNN
      J=1
  660 K1=K1+KSPAN
      K2=K2-KSPAN
      JJ=J
      AK=AA
      BK=BB
      AJ=ZERO
      BJ=ZERO
      K=1
  670 K=K+1
      AK=AT(K)*CK(JJ)+AK
      BK=BT(K)*CK(JJ)+BK
      K=K+1
      AJ=AT(K)*SK(JJ)+AJ
      BJ=BT(K)*SK(JJ)+BJ
      JJ=JJ+J
      IF(JJ .GT. JF) JJ=JJ-JF
      IF(K .LT. JF) GO TO 670
      K=JF-J
      A(K1)=AK-BJ
      B(K1)=BK+AJ
      A(K2)=AK+BJ
      B(K2)=BK-AJ
      J=J+1
      IF(J .LT. K) GO TO 660
      KK=KK+KSPNN
      IF(KK .LE. NN) GO TO 640
      KK=KK-NN
      IF(KK .LE. KSPAN) GO TO 640
C     MULTIPLY BY ROTATION FACTOR (EXCEPT FOR FACTORS OF 2 AND 4)
  700 IF(I .EQ. M) GO TO 800
      KK=JC+1
  710 C2=ONE-CD
      S1=SD
  720 C1=C2
      S2=S1
      KK=KK+KSPAN
  730 AK=A(KK)
      A(KK)=C2*AK-S2*B(KK)
      B(KK)=S2*AK+C2*B(KK)
      KK=KK+KSPNN
      IF(KK .LE. NT) GO TO 730
      AK=S1*S2
      S2=S1*C2+C1*S2
      C2=C1*C2-AK
      KK=KK-NT+KSPAN
      IF(KK .LE. KSPNN) GO TO 730
      C2=C1-(CD*C1+SD*S1)
      S1=S1+(SD*C1-CD*S1)
C     THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C     ERROR. IF ROUNDED ARITHMETIC IS USED, THEY MAY
C     BE DELETED.
C     C1=ZERO/(C2**2+S1**2)+ZERO
C     S1=C1*S1
C     C2=C1*C2
      KK=KK-KSPNN+JC
      IF(KK .LE. KSPAN) GO TO 720
      KK=KK-KSPAN+JC+INC
      IF(KK .LE. JC+JC) GO TO 710
      GO TO 100
C     PERMUTE THE RESULTS TO NORMAL ORDER---DONE IN TWO STAGES
C     PERMUTATION FOR SQUARE FACTORS OF N
  800 NP(1)=KS
      IF(KT .EQ. 0) GO TO 890
      K=KT+KT+1
      IF(M .LT. K) K=K-1
      J=1
      NP(K+1)=JC
  810 NP(J+1)=NP(J)/NFAC(J)
      NP(K)=NP(K+1)*NFAC(J)
      J=J+1
      K=K-1
      IF(J .LT. K) GO TO 810
      K3=NP(K+1)
      KSPAN=NP(2)
      KK=JC+1
      K2=KSPAN+1
      J=1
      IF(N .NE. NTOT) GO TO 850
C     PERMUTATION FOR SINGLE-VARIATE TRANSFORM (OPTIONAL CODE)
  820 AK=A(KK)
      A(KK)=A(K2)
      A(K2)=AK
      BK=B(KK)
      B(KK)=B(K2)
      B(K2)=BK
      KK=KK+INC
      K2=KSPAN+K2
      IF(K2 .LT. KS) GO TO 820
  830 K2=K2-NP(J)
      J=J+1
      K2=NP(J+1)+K2
      IF(K2 .GT. NP(J)) GO TO 830
      J=1
  840 IF(KK .LT. K2) GO TO 820
      KK=KK+INC
      K2=KSPAN+K2
      IF(K2 .LT. KS) GO TO 840
      IF(KK .LT. KS) GO TO 830
      JC=K3
      GO TO 890
C     PERMUTATION FOR MULTIVARIATE TRANSFORM
  850 K=KK+JC
  860 AK=A(KK)
      A(KK)=A(K2)
      A(K2)=AK
      BK=B(KK)
      B(KK)=B(K2)
      B(K2)=BK
      KK=KK+INC
      K2=K2+INC
      IF(KK .LT. K) GO TO 860
      KK=KK+KS-JC
      K2=K2+KS-JC
      IF(KK .LT. NT) GO TO 850
      K2=K2-NT+KSPAN
      KK=KK-NT+JC
      IF(K2 .LT. KS) GO TO 850
  870 K2=K2-NP(J)
      J=J+1
      K2=NP(J+1)+K2
      IF(K2 .GT. NP(J)) GO TO 870
      J=1
  880 IF(KK .LT. K2) GO TO 850
      KK=KK+JC
      K2=KSPAN+K2
      IF(K2 .LT. KS) GO TO 880
      IF(KK .LT. KS) GO TO 870
      JC=K3
  890 IF(2*KT+1 .GE. M) RETURN
      KSPNN=NP(KT+1)
C     PERMUTATION FOR SQUARE-FREE FACTORS OF N
      J=M-KT
      NFAC(J+1)=1
  900 NFAC(J)=NFAC(J)*NFAC(J+1)
      J=J-1
      IF(J .NE. KT) GO TO 900
      KT=KT+1
      NN=NFAC(KT)-1
      IF(NN .GT. MAXP) GO TO 998
      JJ=0
      J=0
      GO TO 906
  902 JJ=JJ-K2
      K2=KK
      K=K+1
      KK=NFAC(K)
  904 JJ=KK+JJ
      IF(JJ .GE. K2) GO TO 902
      NP(J)=JJ
  906 K2=NFAC(KT)
      K=KT+1
      KK=NFAC(K)
      J=J+1
      IF(J .LE. NN) GO TO 904
C     DETERMINE THE PERMUTATION CYCLES OF LENGTH GREATER THAN 1
      J=0
      GO TO 914
  910 K=KK
      KK=NP(K)
      NP(K)=-KK
      IF(KK .NE. J) GO TO 910
      K3=KK
  914 J=J+1
      KK=NP(J)
      IF(KK .LT. 0) GO TO 914
      IF(KK .NE. J) GO TO 910
      NP(J)=-J
      IF(J .NE. NN) GO TO 914
      MAXF=INC*MAXF
C     REORDER A AND B, FOLLOWING THE PERMUTATION CYCLES
      GO TO 950
  924 J=J-1
      IF(NP(J) .LT. 0) GO TO 924
      JJ=JC
  926 KSPAN=JJ
      IF(JJ .GT. MAXF) KSPAN=MAXF
      JJ=JJ-KSPAN
      K=NP(J)
*     KK=JC*K+II+JJ
      KK=JC*K+I+JJ
      K1=KK+KSPAN
      K2=0
  928 K2=K2+1
      AT(K2)=A(K1)
      BT(K2)=B(K1)
      K1=K1-INC
      IF(K1 .NE. KK) GO TO 928
  932 K1=KK+KSPAN
      K2=K1-JC*(K+NP(K))
      K=-NP(K)
  936 A(K1)=A(K2)
      B(K1)=B(K2)
      K1=K1-INC
      K2=K2-INC
      IF(K1 .NE. KK) GO TO 936
      KK=K2
      IF(K .NE. J) GO TO 932
      K1=KK+KSPAN
      K2=0
  940 K2=K2+1
      A(K1)=AT(K2)
      B(K1)=BT(K2)
      K1=K1-INC
      IF(K1 .NE. KK) GO TO 940
      IF(JJ .NE. 0) GO TO 926
      IF(J .NE. 1) GO TO 924
  950 J=K3+1
      NT=NT-KSPNN
*     II=NT-INC+1
      I=NT-INC+1
      IF(NT .GE. 0) GO TO 924
      RETURN
C     ERROR FINISH, INSUFFICIENT ARRAY STORAGE
  998 ISN=0
      PRINT 999
      STOP
  999 FORMAT(44H0ARRAY BOUNDS EXCEEDED WITHIN SUBROUTINE CFT)
      END


      SUBROUTINE NFFT( N )

C CHANGES N INTO THE NEXT INTEGER ALLOWED BY THE FFT ROUTINE
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


