      SUBROUTINE CHKGMX (K,BG,MESH,G2MAX)
C
C  Modules
C
      use precision
#ifdef MPI
      use mpi_siesta
#endif

      IMPLICIT DOUBLE PRECISION (A-H,K,O-Z)

      integer :: i, j, i1, i2, i3

      PARAMETER (ZERO=0.D0,HALF=.5D0,TINY=1.D-8,BIG=1.D20)
      DOUBLE PRECISION K(3),BG(3,3),BM(3,3),G(3)
      INTEGER MESH(3),Node
#ifdef MPI
      integer MPIerror
#endif
C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif
      DO I=1,3
        DO J=1,3
          BM(J,I)=BG(J,I)*MESH(I)
        ENDDO
      ENDDO
      CALL MINVEC (BM,BM)
      GMAX=BIG
      DO I3=-1,1
        DO I2=-1,1
          DO I1=-1,1
            IF (I1.EQ.0.AND.I2.EQ.0.AND.I3.EQ.0) GO TO 20
              G(1)=BM(1,1)*I1+BM(1,2)*I2+BM(1,3)*I3
              G(2)=BM(2,1)*I1+BM(2,2)*I2+BM(2,3)*I3
              G(3)=BM(3,1)*I1+BM(3,2)*I2+BM(3,3)*I3
              GMOD=SQRT(DOT(G,G,3))
              R=HALF*GMOD-DOT(K,G,3)/GMOD
              GMAX=MIN(GMAX,R)
  20        CONTINUE
          ENDDO
        ENDDO
      ENDDO
      IF (GMAX.LT.ZERO) THEN
        if (Node.eq.0) then
          WRITE (6,*) 'CHKGMX: K NOT IN FIRST BZ'
        endif
        STOP
      END IF
      G2MSH=GMAX*GMAX-TINY
      IF (G2MSH.LT.G2MAX) THEN
        if (Node.eq.0) then
*        WRITE(6,*) 'CHKGMX: MESH TOO SPARSE. GMAX HAS BEEN REDUCED'
*        WRITE(6,*) 'CHKGMX: OLD G2MAX =',G2MAX
*        WRITE(6,*) 'CHKGMX: NEW G2MAX =',G2MSH
        endif
        G2MAX=G2MSH
      ENDIF
      END
