C $Id: eispack.f,v 1.1 1999/02/26 10:28:10 wdpgaara Exp $

      SUBROUTINE EISCH1(NM,N,AR,AI,WR,ZR,ZI,IERR,WORK)
C
C     ALL EIGENVALUES AND CORRESPONDING EIGENVECTORS OF A COMPLEX
C     HERMITIAN MATRIX
C
*        ADAPTED TO DOUBLE PRECISION BY J.SOLER 30/10/89
*        FOR SINGLE PRECISION SIMPLY COMMENT OUT NEXT LINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*        PARAMETERS INTRODUCED BY J.SOLER
      PARAMETER (ZERO=0.D0,ONE=1.D0)
      DIMENSION AR(NM,NM),AI(NM,NM),WR(N),ZR(NM,NM),ZI(NM,NM),WORK(*)
      CALL HTRIDI(NM,N,AR,AI,WR,ZI,ZI,WORK)
      DO 100 I=1,N
      DO 50 J=1,N
   50 ZR(I,J)=ZERO
  100 ZR(I,I)=ONE
      CALL TQL2C(NM,N,WR,ZI,ZR,IERR)
      IF(IERR.NE.0) RETURN
      CALL HTRIBK(NM,N,AR,AI,WORK,N,ZR,ZI)
      RETURN
      END

      SUBROUTINE HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
*        ADAPTED TO DOUBLE PRECISION BY J.SOLER 30/10/89
*        FOR SINGLE PRECISION SIMPLY COMMENT OUT NEXT LINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*        PARAMETERS INTRODUCED BY J.SOLER
      PARAMETER (ZERO=0.D0,ONE=1.D0)
      DIMENSION AR(NM,N),AI(NM,N),D(N),E(N),E2(N),TAU(2,N)
*        INTEGER I,J,K,L,N,II,NM,JP1
*        REAL F,FI,G,GI,H,HH,SI,SCALE
      TAU(1,N) = ONE
      TAU(2,N) = ZERO
      DO 100 I = 1, N
  100 D(I) = AR(I,I)
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 130
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(AR(I,K)) + ABS(AI(I,K))
         IF (SCALE .NE. ZERO) GO TO 140
         TAU(1,L) = ONE
         TAU(2,L) = ZERO
  130    E(I) = ZERO
         E2(I) = ZERO
         GO TO 290
  140    DO 150 K = 1, L
            AR(I,K) = AR(I,K) / SCALE
            AI(I,K) = AI(I,K) / SCALE
            H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
  150    CONTINUE
         E2(I) = SCALE * SCALE * H
         G = SQRT(H)
         E(I) = SCALE * G
*           NEXT LINE CHANGED BY J.SOLER.
*        F = CABS(CMPLX(AR(I,L),AI(I,L)))
         F = SQRT(AR(I,L)**2+AI(I,L)**2)
         IF (F .EQ. ZERO) GO TO 160
         TAU(1,L) = (AI(I,L) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
         SI = (AR(I,L) * TAU(2,I) + AI(I,L) * TAU(1,I)) / F
         H = H + F * G
         G = ONE + G / F
         AR(I,L) = G * AR(I,L)
         AI(I,L) = G * AI(I,L)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         AR(I,L) = G
  170    F = ZERO
         DO 240 J = 1, L
            G = ZERO
            GI = ZERO
            DO 180 K = 1, J
               G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
               GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
  180       CONTINUE
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
               G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
               GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
  200       CONTINUE
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * AR(I,J) - TAU(2,J) * AI(I,J)
  240    CONTINUE
         HH = F / (H + H)
         DO 260 J = 1, L
            F = AR(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -AI(I,J)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
            DO 260 K = 1, J
               AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K)
     X                           + FI * TAU(2,K) + GI * AI(I,K)
               AI(J,K) = AI(J,K) - F * TAU(2,K) - G * AI(I,K)
     X                           - FI * E(K) - GI * AR(I,K)
  260    CONTINUE
  270    DO 280 K = 1, L
            AR(I,K) = SCALE * AR(I,K)
            AI(I,K) = SCALE * AI(I,K)
  280    CONTINUE
         TAU(2,L) = -SI
  290    HH = D(I)
         D(I) = AR(I,I)
         AR(I,I) = HH
         AI(I,I) = SCALE * SCALE * H
  300 CONTINUE
      RETURN
      END

      SUBROUTINE TQL2C(NM,N,D,E,Z,IERR)
*        ADAPTED TO DOUBLE PRECISION BY J.SOLER 30/10/89
*        FOR SINGLE PRECISION SIMPLY COMMENT OUT NEXT LINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*        PARAMETERS INTRODUCED BY J.SOLER
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,EPMACH=2.D0**(-23))
      DIMENSION D(N),E(N),Z(NM,N)
*        INTEGER I,J,K,L,M,N,II,NM,MML,IERR
*        VARIABLE MACHEP CHANGED TO PARAMETER EPMACH. J.SOLER.
*        REAL B,C,F,G,H,P,R,S,MACHEP
*        MACHEP=TWO**(-23)
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      F = ZERO
      B = ZERO
      E(N) = ZERO
      DO 240 L = 1, N
         J = 0
*        H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         H = EPMACH * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H
         DO 110 M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
  110    CONTINUE
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
         P = (D(L+1) - D(L)) / (TWO * E(L))
         R = SQRT(P*P+ONE)
         H = D(L) - E(L) / (P + SIGN(R,P))
         DO 140 I = L, N
  140    D(I) = D(I) - H
         F = F + H
         P = D(M)
         C = ONE
         S = ZERO
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+ONE)
            E(I+1) = S * P * R
            S = C / R
            C = ONE / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+ONE)
            E(I+1) = S * E(I) * R
            S = ONE / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
  200    CONTINUE
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
      END

      SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
*        ADAPTED TO DOUBLE PRECISION BY J.SOLER 30/10/89
*        FOR SINGLE PRECISION SIMPLY COMMENT OUT NEXT LINE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*        PARAMETER INTRODUCED BY J.SOLER
      PARAMETER (ZERO=0.D0)
      DIMENSION AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
*        REAL H,S,SI
*        INTEGER I,J,K,L,M,N,NM
      DO 50 K = 1, N
         DO 50 J = 1, M
            ZI(K,J) = (- ZR(K,J)) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
      IF (N .EQ. 1) GO TO 200
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. ZERO) GO TO 140
         DO 130 J = 1, M
            S = ZERO
            SI = ZERO
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
  110       CONTINUE
            S = S / H
            SI = SI / H
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
  200 RETURN
      END


