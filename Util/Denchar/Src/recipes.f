CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    This file contains routines taken or adapted from 
C    'Numerical Recipes, The Art of Scientific Computing' by
C    W.H. Press, S.A. Teukolsky, W.T. Veterling and B.P. Flannery,
C    Cambridge U.P. 1987 and 1992.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C * The routines contained in this file are:
C     SUBROUTINE RATINT
C     SUBROUTINE SPLINE
C     SUBROUTINE SPLINT
C     SUBROUTINE POLINT
C     SUBROUTINE splin
C     SUBROUTINE splinu
C     SUBROUTINE FOUR1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




      SUBROUTINE RATINT(XA,YA,N,X,Y,DY) 
C*****************************************************
C Rational  interpolation 
C Adapted from the Numerical Recipes, 
C Modified for double precision and combined with 
C polinomic interpolation by D. Sanchez-Portal, 1996
C***************************************************** 

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (TINY=1.D-15)
      DIMENSION XA(N),YA(N)
      double precision, dimension(n) :: c, d

      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=DABS(X-XA(I))
        IF (H.LT.TINY)THEN
          Y=YA(I)
          DY=0.0D0
          goto 999
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          IF(DD.EQ.0.0D0)GOTO 100
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE

C*** AS RATIONAL INTERPOLATION DOES NOT CONVERGE,****************** 
C*************WE TRY WITH A POLYNOMIAL ONE*************************


100   CALL POLINT(XA,YA,N,X,Y,DY)

C Exit point
  999 continue

      RETURN
      END





      SUBROUTINE SPLINE(DELT,Y,N,YP1,YPN,Y2) 
C****************************************************** 
C Cubic Spline Interpolation.
C Adapted from Numerical Recipes routines for an uniform 
C grid.
C D. Sanchez-Portal, Oct. 1996.
C*****************************************************

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(N),Y2(N)

      double precision U(N)  ! automatic array
    
      IF (YP1.GT..99D30) THEN
        Y2(1)=0.0D0
        U(1)=0.0D0
      ELSE
        Y2(1)=-0.5D0
        U(1)=(3.0D0/DELT)*((Y(2)-Y(1))/DELT-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=0.5D0
        P=SIG*Y2(I-1)+2.0D0
        Y2(I)=(SIG-1.0D0)/P
        U(I)=(3.0D0*( Y(I+1)+Y(I-1)-2.0D0*Y(I) )/(DELT*DELT)
     *      -SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99D30) THEN
        QN=0.0D0
        UN=0.0D0
      ELSE
        QN=0.5D0
        UN=(3.0D0/DELT)*(YPN-(Y(N)-Y(N-1))/DELT)
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.D0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END






      SUBROUTINE SPLINT(DELT,YA,Y2A,N,X,Y,DYDX) 
C******************************************************
C Cubic Spline Interpolation.
C Adapted from Numerical Recipes routines for an uniform
C grid.
C D. Sanchez-Portal, Oct. 1996.
C*****************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      DIMENSION YA(N),Y2A(N)
      
      NLO=INT(X/DELT)+1
      NHI=NLO+1

 
      IF (DELT.EQ.0.D0) PAUSE 'Bad DELT input.'

      A=NHI-X/DELT-1
      B=1.0D0-A
      Y=A*YA(NLO)+B*YA(NHI)+
     *      ((A**3-A)*Y2A(NLO)+(B**3-B)*Y2A(NHI))*(DELT**2)/6.D0

      DYDX=(YA(NHI)-YA(NLO))/DELT +
     * (-((3*(A**2)-1.D0)*Y2A(NLO))+
     * (3*(B**2)-1.D0)*Y2A(NHI))*DELT/6.D0
      RETURN
      END







      SUBROUTINE POLINT(XA,YA,N,X,Y,DY) 
C*****************************************************
C Polinomic interpolation
C Adapted from Numerical Recipes for double precision,
C D. Sanchez-Portal, Oct. 1996
C*****************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(N),YA(N)

      double precision, dimension(n) :: c, d

      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.0D0)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE

      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software !E#.



      SUBROUTINE splin(x,y,n,yp1,ypn,y2)
C This is routine spline of Num. Recipes with the name changed.
C J.M.Soler, April'96.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION x(n),y(n),y2(n)

      double precision, dimension(n) ::  u  ! automatic array

      if (yp1.gt..99d30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.0d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      DO 20 I=1,N
C        WRITE(*,*)'x(',I,') = ',x(I)
C        WRITE(*,*)'y(',I,') = ',y(I)
C        WRITE(*,*)'y2(',I,') = ',y2(I)
20    CONTINUE

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software !E#.




      SUBROUTINE splinu(xmin,xmax,ya,y2a,n,x,y,dy)
C Adapted from the routine splint of Num. Recipes for a uniform grid.
C J.M.Soler, April'96.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION y2a(n),ya(n)
      h=(xmax-xmin)/(n-1)
      klo=(x-xmin)/h+1
      klo=max(1,klo)
      klo=min(klo,n-1)
      khi=klo+1
      hy=(ya(khi)-ya(klo))/h
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xmin+(khi-1)*h-x)/h
      b=(x-xmin-(klo-1)*h)/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      dy=(h*y2a(klo)/6.d0)*((-3.d0)*a**2+1.d0)+
     *   (h*y2a(khi)/6.d0)*(3.d0*b**2-1.d0)+hy

      return
      END


      SUBROUTINE FOUR1(DATA,NN,ISIGN)
C Converted to double precision from same routine in
C "Numerical Recipes", W.Press et al, Cambridge U.P., 1st ed.
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=(-2.D0)*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END


