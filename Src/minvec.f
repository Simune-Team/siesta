C $Id: minvec.f,v 1.6 1999/03/04 16:21:09 jose Exp $

      SUBROUTINE MINVEC(B0,BMIN)

C *******************************************************************
C  FINDS THE LATTICE BASIS OF MINIMUM LENGTH, I.E. SUCH TAHT ANY 
C  OTHER BASIS (NOT EQUIVALENT BY SYMMETRY) HAS ONE VECTOR LONGER.
C  WRITTEN BY J.MORENO AND J.SOLER. AUGUST 1989 AND OCTOBER 1997.
C ********* INPUT ***************************************************
C REAL*8 B0 : Cell vectors B0(xyj,vector)
C ********* OUTPUT **************************************************
C REAL*8 BMIN : Minimum cell vectors B0(xyj,vector)
C *******************************************************************

      IMPLICIT         NONE
      DOUBLE PRECISION B0(3,3),BMIN(3,3),DOT,VOLCEL
      EXTERNAL         DOT,ORDIX,ORDER,RECLAT,VOLCEL,VOLNEW

      INTEGER          I,I1,I2,I3,IAUX(3),ITER,J,NITER
      DOUBLE PRECISION AUX(3,3),B(3,3),B2(3),BNEW(3),BNEW2,
     .                 C(3,3),EPS,VNEW,V0

      PARAMETER (EPS=1.D-8,NITER=100)

      V0=ABS(VOLCEL(B0))
      IF (V0.LT.EPS) THEN
        WRITE(6,*) 'MINVEC: BASIS VECTORS ARE LINEARLY DEPENDENT'
        STOP
      ENDIF
      DO 20 I=1,3
        DO 10 J=1,3
          B(J,I)=B0(J,I)
  10    CONTINUE
        B2(I)=DOT(B(1,I),B(1,I),3)
  20  CONTINUE

      DO 50 ITER=1,NITER
        CALL ORDIX(B2,1,3,IAUX)
        CALL ORDER(B2,1,3,IAUX)
        CALL ORDER(B ,3,3,IAUX)
        DO 40 I1=0,1
        DO 40 I2=-1,1
        DO 40 I3=-1,1
          IF (I1.EQ.0.AND.I2.NE.1) GO TO 40
          IF (I2.EQ.0.AND.I3.EQ.0) GO TO 40
          BNEW(1)=B(1,1)*I1+B(1,2)*I2+B(1,3)*I3
          BNEW(2)=B(2,1)*I1+B(2,2)*I2+B(2,3)*I3
          BNEW(3)=B(3,1)*I1+B(3,2)*I2+B(3,3)*I3
          BNEW2=DOT(BNEW,BNEW,3)
          DO 30 I=3,1,-1
            IF (BNEW2+EPS.GE.B2(I)) GO TO 40
              CALL VOLNEW(B,BNEW,I,VNEW)
              IF (ABS((VNEW-V0)/V0).LT.EPS) THEN
                B(1,I)=BNEW(1)
                B(2,I)=BNEW(2)
                B(3,I)=BNEW(3)
                B2(I)=BNEW2
                GO TO 50
              END IF
  30      CONTINUE
  40    CONTINUE
        GOTO 55
  50  CONTINUE
        WRITE(6,*) 'MINVEC: ERROR: Iteration has not converged'
        STOP 'MINVEC: ERROR: Iteration has not converged'
  55  CONTINUE

      IF (VOLCEL(B).LT.0.D0) THEN
        B(1,3)=-B(1,3)
        B(2,3)=-B(2,3)
        B(3,3)=-B(3,3)
      ENDIF
      CALL RECLAT(B0,AUX,0)
      DO 60 I=1,3
      DO 60 J=1,3
        C(J,I)=NINT(DOT(AUX(1,J),B(1,I),3))
  60  CONTINUE
      DO 70 I=1,3
      DO 70 J=1,3
        B(J,I)=B0(J,1)*C(1,I)+B0(J,2)*C(2,I)+B0(J,3)*C(3,I)
  70  CONTINUE
      DO 80 I=1,3
      DO 80 J=1,3
        BMIN(J,I)=B(J,I)
  80  CONTINUE
      END


      SUBROUTINE VOLNEW(A,ANEW,INEW,VOL)
      IMPLICIT         NONE
      INTEGER          I,INEW,J
      DOUBLE PRECISION A(3,3),ANEW(3),AUX(3,3),VOL,VOLCEL
      DO 10 I=1,3
      DO 10 J=1,3
        AUX(J,I)=A(J,I)
  10  CONTINUE
      DO 20 J=1,3
        AUX(J,INEW)=ANEW(J)
  20  CONTINUE
      VOL=ABS(VOLCEL(AUX))
      END
