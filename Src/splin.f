      SUBROUTINE splin(x,y,n,yp1,ypn,y2)
C This is routine spline of Num. Recipes with the name changed.
C J.M.Soler, April'96.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=2000)
      DIMENSION x(n),y(n),y2(n),u(NMAX)
      if (n.gt.NMAX) stop 'splin: Dimension NMAX too small'
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
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
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xmin+(khi-1)*h-x)/h
      b=(x-xmin-(klo-1)*h)/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      dy=(h*y2a(klo)/6.)*((-3.)*a**2+1.)+
     *   (h*y2a(khi)/6.)*(3.*b**2-1.)+hy

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software !E#.

