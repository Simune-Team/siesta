      module m_recipes
!**********************************************************************
!    This file contains routines adapted from 'Numerical Recipes, 
!    The Art of Scientific Computing' by W.H. Press, S.A. Teukolsky, 
!    W.T. Veterling and B.P. Flannery, Cambridge U.P. 1987 and 1992.
!**********************************************************************

! Used module procedures and types
      use sys,       only: die  ! Termination routine
      use precision, only: dp   ! Double precision real kind

      implicit none

! Public procedures provided by this module
      public :: four1   ! 1-D fast Fourier transform
      public :: tqli    ! With tred2, diagonalizes a real matrix
      public :: tred2   ! Reduction of a real matrix to tridiagonal form

      private ! Nothing is declared public below this point

      contains

      SUBROUTINE FOUR1(DATA,NN,ISIGN)
!**********************************************************************
! Discrete Fourier transform. Modified and converted to double 
! precision from same routine in Numerical Recipes.
!**********************************************************************
! Input:
!   complex*16 DATA(NN) : Function to be Fourier transformed
!   integer    NN       : Number of points. Must be a power of 2
!   integer    ISIGN    : ISIG=+1/-1 => Direct/inverse transform
! Output:
!   complex*16 DATA(NN) : The direct Fourier transform (ISIG=+1), or
!                         NN times the inverse Fourier transf (ISIG=-1)
!**********************************************************************
      IMPLICIT NONE
      INTEGER  :: NN, ISIGN
      REAL(dp) :: DATA(2*NN)

      INTEGER  :: I, ISTEP, J, M, MMAX, N
      REAL(dp) :: TEMPI, TEMPR, THETA, WI, WPI, WPR, WR, WTEMP
      REAL(dp), PARAMETER :: TWOPI=6.28318530717959D0,
     .  HALF=0.5D0, ONE=1.D0, TWO=2.D0, ZERO=0.D0

      N=2*NN
      J=1
      DO I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
        DO ! until following condition is met
          IF ((M.LT.2).OR.(J.LE.M)) EXIT
          J=J-M
          M=M/2
        END DO
        J=J+M
      END DO ! I
      MMAX=2
      DO ! until following condition is met
        IF (N.LE.MMAX) EXIT
        ISTEP=2*MMAX
        THETA=TWOPI/(ISIGN*MMAX)
        WPR=(-TWO)*SIN(HALF*THETA)**2
        WPI=SIN(THETA)
        WR=ONE
        WI=ZERO
        DO M=1,MMAX,2
          DO I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
          END DO ! I
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
        END DO ! M
        MMAX=ISTEP
      END DO ! until (N.LE.MMAX)

      END SUBROUTINE FOUR1


      SUBROUTINE TQLI(D,E,N,NP,Z)

! IN COMBINATION WITH TRED2 FINDS EIGENVALUES AND EIGENVECTORS OF
! A REAL SYMMETRIC MATRIX. REF: W.H.PRESS ET AL. NUMERICAL RECIPES.

      implicit none

      integer, intent(in)   :: N        ! True matrix dimension
      integer, intent(in)   :: NP       ! Physical size of arrays
      real(dp),intent(inout):: D(NP)    ! Vector prepared by tred2
      real(dp),intent(inout):: E(NP)    ! Input: vector prepared by tred2
                                        ! Output: eigenvalues
      real(dp),intent(out)  :: Z(NP,NP) ! Eigenvectors, in columns
      
      real(dp) :: zero, one, two
      PARAMETER (ZERO=0.0_dp,ONE=1.0_dp,TWO=2.0_dp)

      integer  :: iter, i, k, l, m
      real(dp) :: dd, g, r, s, c, p, f, b

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
      END subroutine tqli


      SUBROUTINE TRED2(A,N,NP,D,E)

!HOUSEHOLDER REDUCTION OF A REAL SYMMETRIC MATRIX INTO TRIDIAGONAL FORM
!REF: W.H.PRESS ET AL. NUMERICAL RECIPES. CAMBRIDGE U.P.

      implicit none

      integer, intent(in)   :: N        ! True matrix dimension
      integer, intent(in)   :: NP       ! Physical size of arrays
      real(dp),intent(inout):: A(NP,NP) ! Real symmetric matrix (overwritten)
      real(dp),intent(out)  :: D(NP)    ! Vector to be used by tqli
      real(dp),intent(out)  :: E(NP)    ! Vector to be used by tqli

      real(dp) ::  zero, one
      PARAMETER (ZERO=0.0_dp,ONE=1.0_dp)

      integer  :: l, i, k, j
      real(dp) :: f, g, h, hh, scale

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
      END subroutine tred2

      end module m_recipes
