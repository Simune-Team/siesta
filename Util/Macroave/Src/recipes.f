!**********************************************************************
!    This file contains routines adapted from 'Numerical Recipes, 
!    The Art of Scientific Computing' by W.H. Press, S.A. Teukolsky, 
!    W.T. Veterling and B.P. Flannery, Cambridge U.P. 1987 and 1992.
!**********************************************************************
! The routines contained in this file are:
!     SUBROUTINE FOUR1
!**********************************************************************


      SUBROUTINE FOUR1(DATA,NN,ISIGN)
!**********************************************************************
! Discrete Fourier transform. Modified and converted to double 
! precision from same routine in Numerical Recipes.
!**********************************************************************
! Input:
!   real*8  DATA(NN) : Function to be Fourier transformed
!   integer NN       : Number of points. Must be a power of 2
!   integer ISIGN    : ISIG=+1/-1 => Direct/inverse transform
! Output:
!   real*8  DATA(NN) : Fourier transformed function
!**********************************************************************
      IMPLICIT NONE
      INTEGER          :: NN, ISIGN
      DOUBLE PRECISION :: DATA(NN)

      INTEGER          :: I, ISTEP, J, M, MMAX, N
      DOUBLE PRECISION :: TEMPI, TEMPR, THETA, WI, WPI, WPR, WR, WTEMP
      DOUBLE PRECISION, PARAMETER :: TWOPI=6.28318530717959D0,
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

