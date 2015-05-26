      module m_recipes
!**********************************************************************
!    This file contains routines adapted from 'Numerical Recipes, 
!    The Art of Scientific Computing' by W.H. Press, S.A. Teukolsky, 
!    W.T. Veterling and B.P. Flannery, Cambridge U.P. 1987 and 1992.
!**********************************************************************

! Used module procedures and types
      use precision, only: dp   ! Double precision real kind

      implicit none

! Public procedures provided by this module
      public :: four1   ! 1-D fast Fourier transform

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


      end module m_recipes
