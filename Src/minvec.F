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

C
C  Modules
C
      use precision
#ifdef MPI
      use mpi
#endif

      IMPLICIT         NONE

#ifdef MPI
      integer          MPIerror
#endif
      DOUBLE PRECISION B0(3,3),BMIN(3,3),DOT,VOLCEL
      EXTERNAL         DOT,ORDIX,ORDER,RECLAT,VOLCEL,VOLNEW

      INTEGER          I,I1,I2,I3,ITER,J,NITER,Node,IAUX(3)
      DOUBLE PRECISION AUX(3,3),B(3,3),B2(3),BNEW(3),BNEW2,
     .                 C(3,3),EPS,VNEW,V0

      PARAMETER (EPS=1.D-8,NITER=100)

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif
      V0=ABS(VOLCEL(B0))
      IF (V0.LT.EPS) THEN
        if (Node.eq.0) then
          WRITE(6,*) 'MINVEC: BASIS VECTORS ARE LINEARLY DEPENDENT'
        endif
        STOP
      ENDIF
      DO I=1,3
        DO J=1,3
          B(J,I)=B0(J,I)
        ENDDO
        B2(I)=DOT(B(1,I),B(1,I),3)
      ENDDO

      DO 50 ITER=1,NITER
        CALL ORDIX(B2,1,3,IAUX)
        CALL ORDER(B2,1,3,IAUX)
        CALL ORDER(B ,3,3,IAUX)
        DO I1=0,1
          DO I2=-1,1
            DO I3=-1,1
              IF (I1.EQ.0.AND.I2.NE.1) GO TO 40
              IF (I2.EQ.0.AND.I3.EQ.0) GO TO 40
              BNEW(1)=B(1,1)*I1+B(1,2)*I2+B(1,3)*I3
              BNEW(2)=B(2,1)*I1+B(2,2)*I2+B(2,3)*I3
              BNEW(3)=B(3,1)*I1+B(3,2)*I2+B(3,3)*I3
              BNEW2=DOT(BNEW,BNEW,3)
              DO I=3,1,-1
                IF (BNEW2+EPS.GE.B2(I)) GO TO 40
                CALL VOLNEW(B,BNEW,I,VNEW)
                IF (ABS((VNEW-V0)/V0).LT.EPS) THEN
                  B(1,I)=BNEW(1)
                  B(2,I)=BNEW(2)
                  B(3,I)=BNEW(3)
                  B2(I)=BNEW2
                  GO TO 50
                END IF
              ENDDO
  40          CONTINUE
            ENDDO
          ENDDO
        ENDDO
        GOTO 55
  50  CONTINUE
      if (Node.eq.0) then
        WRITE(6,*) 'MINVEC: ERROR: Iteration has not converged'
        STOP 'MINVEC: ERROR: Iteration has not converged'
      else
        STOP
      endif
  55  CONTINUE

      IF (VOLCEL(B).LT.0.D0) THEN
        B(1,3)=-B(1,3)
        B(2,3)=-B(2,3)
        B(3,3)=-B(3,3)
      ENDIF
      CALL RECLAT(B0,AUX,0)
      DO I=1,3
        DO J=1,3
          C(J,I)=NINT(DOT(AUX(1,J),B(1,I),3))
        ENDDO
      ENDDO
      DO I=1,3
        DO J=1,3
          B(J,I)=B0(J,1)*C(1,I)+B0(J,2)*C(2,I)+B0(J,3)*C(3,I)
        ENDDO
      ENDDO
      DO I=1,3
        DO J=1,3
          BMIN(J,I)=B(J,I)
        ENDDO
      ENDDO
      END


      SUBROUTINE VOLNEW(A,ANEW,INEW,VOL)
      IMPLICIT         NONE
      INTEGER          I,INEW,J
      DOUBLE PRECISION A(3,3),ANEW(3),AUX(3,3),VOL,VOLCEL
      AUX(1:3,1:3)=A(1:3,1:3)
      AUX(1:3,INEW)=ANEW(1:3)
      VOL=ABS(VOLCEL(AUX))
      END
