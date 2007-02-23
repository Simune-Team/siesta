!
!  Helper routines for filtering package:
!
!  ordix.f order.f plgndr.f tqli.f tred2.f xlgndr.f
!
!  Note that the "order" routines in module 'sorting' are not used by
!  filter, as it does not use that module. It uses the versions here
!  instead.
!
      SUBROUTINE ORDIX( X, M, N, INDX )
C *******************************************************************
C Makes an index table of array X, with size N and stride M.
C Ref: W.H.Press et al. Numerical Recipes, Cambridge Univ. Press.
C Adapted by J.M.Soler from routine INDEXX of Num. Rec. May'96.
C *************** INPUT *********************************************
C REAL*8  X(M,N)   : Array with the values to be ordered
C INTEGER M, N     : Dimensions of array X
C *************** OUTPUT ********************************************
C INTEGER INDEX(N) : Array which gives the increasing order of X(1,I):
C                    X(1,INDEX(I)) .LE. X(1,INDEX(I+1)) )
C *************** USAGE *********************************************
C Example to order atomic positions X(I,IA), I=1,3, IA=1,NA by
C increasing z coordinate:
C    CALL ORDIX( X(3,1), 3, NA, INDEX )
C    CALL ORDER( X(1,1), 3, NA, INDEX )
C *******************************************************************
      IMPLICIT          NONE
      INTEGER           I, N, INDX(N), INDXT, IR, J, L, M
      DOUBLE PRECISION  X(M,N), Q

      DO 1 J=1,N
         INDX(J)=J
   1  CONTINUE
      IF (N.LE.1) RETURN
      L=N/2+1
      IR=N
   2  CONTINUE
         IF (L.GT.1) THEN
            L=L-1
            INDXT=INDX(L)
            Q=X(1,INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=X(1,INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF (IR.EQ.1) THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
   3     IF (J.LE.IR) THEN
            IF (J.LT.IR) THEN
               IF (X(1,INDX(J)).LT.X(1,INDX(J+1))) J=J+1
            ENDIF
            IF (Q.LT.X(1,INDX(J))) THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GO TO 3
         ENDIF
         INDX(I)=INDXT
      GO TO 2

      END


      SUBROUTINE ORDER( X, M, N, INDEX )
C *******************************************************************
C Orders array X(M,N) according to array INDEX(N), which may be
C   generated by routine ORDIX.
C Written by J.M.Soler. May'96.
C *************** INPUT *********************************************
C INTEGER M, N     : Dimensions of array X
C INTEGER INDEX(N) : Array which gives the desired order
C *************** INPUT AND OUTPUT **********************************
C REAL*8  X(M,N) : Array(s) to be ordered: Xout(I,J) = Xin(I,INDEX(J))
C *******************************************************************
      IMPLICIT          NONE
      INTEGER           I, N, INDEX(N), IORDER, ISTORE, J, M
      DOUBLE PRECISION  X(M,N), XI

      DO 40 J = 1,M
        DO 20 I = 1,N
          XI = X(J,I)
          IORDER = I
   10     CONTINUE
          ISTORE = INDEX(IORDER)
          IF (ISTORE .GT. 0) THEN
            IF (ISTORE .EQ. I) THEN
              X(J,IORDER) = XI
            ELSE
              X(J,IORDER) = X(J,ISTORE)
            ENDIF
            INDEX(IORDER) = -INDEX(IORDER)
            IORDER = ISTORE
            GOTO 10
          ENDIF
   20   CONTINUE
        DO 30 I = 1,N
          INDEX(I) = -INDEX(I)
   30   CONTINUE
   40 CONTINUE
      END
      DOUBLE PRECISION FUNCTION PLGNDR(L,M,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ONE=1.D0,TWO=2.D0)
      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.ONE)STOP ' PLGNDR: bad arguments'
      PMM=ONE
      IF(M.GT.0) THEN
        SOMX2=SQRT((ONE-X)*(ONE+X))
        FACT=ONE
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+TWO
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE TQLI(D,E,N,NP,Z)
C
C  IN COMBINATION WITH TRED2 FINDS EIGENVALUES AND EIGENVECTORS OF
C  A REAL SYMMETRIC MATRIX. REF: W.H.PRESS ET AL. NUMERICAL RECIPES.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0)
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=ZERO
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30) STOP 'tqli: too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(TWO*E(L))
            R=SQRT(G**2+ONE)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=ONE
            C=ONE
            P=ZERO
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+ONE)
                E(I+1)=F*R
                S=ONE/R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+ONE)
                E(I+1)=G*R
                C=ONE/R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+TWO*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=ZERO
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE TRED2(A,N,NP,D,E)
C
C HOUSEHOLDER REDUCTION OF A REAL SYMMETRIC MATRIX INTO TRIDIAGONAL FORM
C REF: W.H.PRESS ET AL. NUMERICAL RECIPES. CAMBRIDGE U.P.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NP,NP),D(NP),E(NP)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=ZERO
          SCALE=ZERO
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.ZERO)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=ZERO
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=ZERO
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=ZERO
      E(1)=ZERO
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.ZERO)THEN
          DO 21 J=1,L
            G=ZERO
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=ONE
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=ZERO
            A(J,I)=ZERO
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END
      SUBROUTINE XLGNDR (X,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),W(N)
      PARAMETER (EPS=3.D-14,PI=3.14159265358979D0)
      M=(N+1)/2
      DO 12 I=1,M
        Z=COS(PI*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        J=M+1-I
        X(J)=Z
        W(J)=2.D0/((1.D0-Z*Z)*PP*PP)
        IF(ABS(Z).LT.EPS) W(J)=.5D0*W(J)
12    CONTINUE
      RETURN
      END