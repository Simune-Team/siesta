C $Id: cdiag_general.f,v 1.2 1999/03/11 19:29:14 ordejon Exp $

      SUBROUTINE CDIAG (H,NH,S,NS,N,E,Z,NZ,AUX)

*  SOLVES THE COMPLEX GENERAL EIGENVALUE PROBLEM   H*Z = E*S*Z
*  WHERE S AND H ARE COMPLEX HERMITIAN MATRICES

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(2*NH*N),S(2*NS*N),Z(2*NZ*N),E(N),AUX(N,5)

*  START TIME COUNT
      CALL TIMER('CDIAG',1)

*  REDUCE GENERAL EIGENVALUE PROBLEM TO A NORMAL ONE
      CALL REDUCC (NH,NS,N,H,S,AUX)

*  NEXT LINE FOR IMSL
C     CALL EIGCH (H,NH,11,E,Z,NZ,AUX(1,2),IER)

*  NEXT LINE FOR SLATEC. AUX 4*N REALS
C     CALL CHIEV (H,NH,N,E,Z,NZ,AUX(1,2),1,INFO)

*  NEXT THREE LINES FOR IBM-ESSL. AUX 4*N REALS
C     CALL PACK (NH,N,H,H,+1)
C     CALL ZHPEV (1,H,E,Z,NZ,N,AUX(1,2),4*N)
C     CALL PACK (NH,N,H,H,-1)

*  NEXT LINES FOR EISPACK. AUX 2*N REALS.
      IF (NH.NE.NZ) THEN
        WRITE (6,*) 'CDIAG: NH AND NZ MUST BE EQUAL. NH,NZ =',NH,NZ
        STOP 2
      ENDIF
      CALL TRANS (2,NH*N,H,H,Z)
      CALL EISCH1 (NH,N,H,H(NH*N+1),E,Z,Z(NZ*N+1),IER,AUX(1,2))
      IF (IER.NE.0) THEN
        WRITE (6,*) 'CDIAG: AN ERROR OCCURRED IN EISCH1. IER =',IER
        STOP 3
      ENDIF
      CALL TRANS (NZ*N,2,Z,Z,H)

*  UNDO THE REDUCTION OF GENERAL EIGENVALUE PROBLEM BY REDUCC.
      CALL REBAKK (NS,NZ,N,S,AUX,Z)

*  STOP TIME COUNT
      CALL TIMER('CDIAG',2)
      END


      SUBROUTINE TRANS (N1,N2,A,B,AUX)

*  TRANSPOSES MATRIX A(N1,N2) TO B(N2,N1). A AND B MAY BE
*  THE SAME PHYSICAL ARRAY. WRITTEN BY J.SOLER

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N1,N2),B(N2,N1),AUX(N1*N2)
      K=0
      DO 20 I=1,N2
        DO 10 J=1,N1
          K=K+1
          AUX(K)=A(J,I)
   10   CONTINUE
   20 CONTINUE
      K=0
      DO 40 I=1,N2
        DO 30 J=1,N1
          K=K+1
          B(I,J)=AUX(K)
   30   CONTINUE
   40 CONTINUE
      END


      SUBROUTINE PACK (NA,N,A,AP,ISN)

*  PACKS (ISN=1) OR UNPACKS (ISN=-1) A COMPLEX HERMITIAN
*  MATRIX TO OR FROM PACKED LOWER MODE.  WRITTEN BY J.SOLER

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(2,NA,N),AP(2,*)
      PARAMETER (ZERO=0.D0)
      IF (ISN.GT.0) THEN
         K=N+1
         DO 20 I=2,N
            DO 10 J=I,N
               AP(1,K)=A(1,J,I)
               AP(2,K)=A(2,J,I)
               K=K+1
  10        CONTINUE
  20     CONTINUE
      ELSE
         K=(N*(N+1))/2
         DO 40 I=N,2,-1
            DO 30 J=N,I,-1
               A(1,J,I)=AP(1,K)
               A(2,J,I)=AP(2,K)
               K=K-1
  30        CONTINUE
  40     CONTINUE
         DO 60 I=1,N
            DO 50 J=1,I-1
               A(1,J,I)= A(1,I,J)
               A(2,J,I)=-A(2,I,J)
  50        CONTINUE
            A(2,I,I)=ZERO
  60     CONTINUE
      ENDIF
      END


      SUBROUTINE REDUCC (NH,NS,N,H,S,DL)

*  GENERALIZATION OF REDUC1 (COMPLEX)
*  CHOLESKY FACTORIZATION OF S ( S=L*TRANSPOSE(L) )
*    USING SCHMIDT' ORTHONORMALIZATION METHOD
*    I)  REPRESENTS A NON-ORTHOGONAL BASIS FUNCTION
*    I>  REPRESENTS A ORTHONORMAL WAVE FUNCTION
*  OBTAINED FROM M.METHFESSEL. SLIGHTLY REWRITTEN BY J.SOLER.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(2,NH,N),S(2,NS,N),DL(N)
      PARAMETER (ZERO=0.D0,ONE=1.D0)

*  OBTAIN TRANSFER MATRIX <IJ) FOR I<J USING CHOLESKY FACTORIZATION
      DO 30 I=1,N
         Y=ONE/DL(I)
         DO 20 J=I,N
            XR=S(1,I,J)
            XI=S(2,I,J)
            DO 10 K=1,I-1
               XR=XR-S(1,J,K)*S(1,I,K)-S(2,J,K)*S(2,I,K)
               XI=XI+S(1,J,K)*S(2,I,K)-S(2,J,K)*S(1,I,K)
  10        CONTINUE
            IF (J.EQ.I) THEN
               IF (XR.LE.ZERO) THEN
                  WRITE (6,*) 'REDUCC: MATRIX S NOT POS. DEFINITE. I=',I
                  STOP
               END IF
               DL(I)=DSQRT(XR)
               Y=ONE/DL(I)
            ELSE
               S(1,J,I)=XR*Y
               S(2,J,I)=XI*Y
            END IF
  20     CONTINUE
  30  CONTINUE

*  OBTAIN MATRIX <IHJ) FOR I.LE.J
      DO 60 I=1,N
         Y=ONE/DL(I)
         DO 50 J=I,N
            XR=H(1,J,I)
            XI=H(2,J,I)
            DO 40 K=1,I-1
               XR=XR-H(1,J,K)*S(1,I,K)+H(2,J,K)*S(2,I,K)
               XI=XI-H(1,J,K)*S(2,I,K)-H(2,J,K)*S(1,I,K)
  40        CONTINUE
            H(1,J,I)=XR*Y
            H(2,J,I)=XI*Y
  50     CONTINUE
  60  CONTINUE

*  OBTAIN MATRIX <IHJ> FOR I.LE.J
      DO 100 I=1,N
         DO 90 J=I,N
            XR=H(1,J,I)
            XI=H(2,J,I)
            DO 70 K=1,I-1
               XR=XR-S(1,J,K)*H(1,I,K)+S(2,J,K)*H(2,I,K)
               XI=XI+S(1,J,K)*H(2,I,K)+S(2,J,K)*H(1,I,K)
  70        CONTINUE
            DO 80 K=I,J-1
               XR=XR-S(1,J,K)*H(1,K,I)-S(2,J,K)*H(2,K,I)
               XI=XI-S(1,J,K)*H(2,K,I)+S(2,J,K)*H(1,K,I)
  80        CONTINUE
            H(1,J,I)=XR/DL(J)
            H(2,J,I)=XI/DL(J)
  90     CONTINUE
 100  CONTINUE

*  OBTAIN MATRIX <IHJ> FOR I.GE.J
      DO 120 I=1,N
         H(2,I,I)=ZERO
         DO 110 J=I+1,N
            H(1,I,J)= H(1,J,I)
            H(2,I,J)=-H(2,J,I)
 110     CONTINUE
 120  CONTINUE
      END


      SUBROUTINE REBAKK (NS,NZ,N,S,DL,Z)

*  BACK TRANSFORMATION OF THE EIGENVECTORS OF THE
*  GENERAL EIGENVALUE PROBLEM:  A*X = LAMBDA*B*X
*  THE FORWARD TRANSFORMATION WAS PERFORMED BY ROUTINE 'REDUCC'

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(2,NS,N),Z(2,NZ,N),DL(N)
      DO 30 I=1,N
         DO 20 J=N,1,-1
            XR=Z(1,J,I)
            XI=Z(2,J,I)
            DO 10 K=J+1,N
               XR=XR-S(1,K,J)*Z(1,K,I)+S(2,K,J)*Z(2,K,I)
               XI=XI-S(1,K,J)*Z(2,K,I)-S(2,K,J)*Z(1,K,I)
  10        CONTINUE
            Z(1,J,I)=XR/DL(J)
            Z(2,J,I)=XI/DL(J)
  20     CONTINUE
  30  CONTINUE
      END




