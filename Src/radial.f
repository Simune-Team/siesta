      module radial

      use precision
      use types

      implicit none

      real(dp)             :: rad_rvals(nrtmax)
      real(dp)             :: rad_fvals(nrtmax)

      private 
      public rad_get, rad_setup_d2, rad_zero

      contains

      subroutine rad_get(func,r,fr,dfdr)
      type(rad_func), intent(in) :: func
      real(dp), intent(in)         :: r
      real(dp), intent(out)        :: fr
      real(dp), intent(out)        :: dfdr

      external splint

      call splint(func%delta,func%f,func%d2,nrtmax,r,fr,dfdr)
      
      end subroutine rad_get
!
!     Set up second derivative in a radial function
!
      subroutine rad_setup_d2(func)
      type(rad_func), intent(inout) :: func

      real(dp) aux(nrtmax)
      real(dp) yp1, ypn
      external spline

      yp1 = huge(1._dp)
      ypn = huge(1._dp)
      call spline(func%delta,func%f,nrtmax,yp1,ypn,func%d2,aux)
      
      end subroutine rad_setup_d2

      subroutine rad_zero(func)
      type(rad_func), intent(inout) :: func

      func%delta  = 0._dp
      func%cutoff = 0._dp
      func%f(:)   = 0._dp
      func%d2(:)  = 0._dp
      
      end subroutine rad_zero
!
!     Fill in a radial function
!
      subroutine rad_put(func,delta,cutoff,vals)
      type(rad_func), intent(inout) :: func
      real(dp), intent(in)         :: delta
      real(dp), intent(in)         :: cutoff
      real(dp), intent(in)         :: vals(:)

      real(dp) aux(nrtmax)
      real(dp) yp1, ypn
      integer i

      external spline

      yp1 = huge(1._dp)
      ypn = huge(1._dp)
      func%delta = delta
      func%cutoff = cutoff
      do i=1,nrtmax-1
         func%f(i) = vals(i)
      enddo

      call spline(delta,func%f,nrtmax,yp1,ypn,func%d2,aux)
      
      end subroutine rad_put
!
!
!
      subroutine rad_getrvals(func)
      type(rad_func), intent(in) :: func

      integer i

      do i=1,nrtmax
         rad_rvals(i) = func%delta *(i-1)
      enddo
      end subroutine rad_getrvals


      SUBROUTINE SPLINT(DELT,YA,Y2A,N,X,Y,DYDX) 
C Cubic Spline Interpolation.
C Adapted from Numerical Recipes for a uniform grid.

      implicit none

      real(dp), intent(in)  :: delt
      real(dp), intent(in)  :: ya(:), y2a(:)
      integer, intent(in)   :: n
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: y
      real(dp), intent(out) :: dydx
      
      integer nlo, nhi
      real(dp) a, b

      NLO=INT(X/DELT)+1
      NHI=NLO+1

      A=NHI-X/DELT-1
      B=1.0_DP-A

      Y=A*YA(NLO)+B*YA(NHI)+
     *      ((A**3-A)*Y2A(NLO)+(B**3-B)*Y2A(NHI))*(DELT**2)/6._DP

      DYDX=(YA(NHI)-YA(NLO))/DELT +
     * (-((3*(A**2)-1._DP)*Y2A(NLO))+
     * (3*(B**2)-1._DP)*Y2A(NHI))*DELT/6._DP

      end subroutine splint

      SUBROUTINE SPLINE(DELT,Y,N,YP1,YPN,Y2,U) 

C Cubic Spline Interpolation.
C Adapted from Numerical Recipes routines for a uniform grid
C D. Sanchez-Portal, Oct. 1996.
! Alberto Garcia, June 2000

      implicit none

      integer, intent(in)    :: n
      real(dp), intent(in)   :: delt
      real(dp), intent(in)   :: yp1, ypn
      real(dp), intent(in)   :: Y(:)
      real(dp), intent(out)  :: Y2(:),U(:)

      integer i, k
      real(dp) sig, p, qn, un

      IF (YP1.eq. huge(1._dp)) THEN
        Y2(1)=0.0_dp
        U(1)=0.0_dp
      ELSE
        Y2(1)=-0.5_dp
        U(1)=(3.0_DP/DELT)*((Y(2)-Y(1))/DELT-YP1)
      ENDIF

      DO I=2,N-1
        SIG=0.5_DP
        P=SIG*Y2(I-1)+2.0_DP
        Y2(I)=(SIG-1.0_DP)/P
        U(I)=(3.0_DP*( Y(I+1)+Y(I-1)-2.0_DP*Y(I) )/(DELT*DELT)
     $      -SIG*U(I-1))/P
      ENDDO

      IF (YPN.eq.huge(1._dp)) THEN
        QN=0.0_DP
        UN=0.0_DP
      ELSE
        QN=0.5_DP
        UN=(3.0_DP/DELT)*(YPN-(Y(N)-Y(N-1))/DELT)
      ENDIF

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1._DP)
      DO K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      enddo

      END subroutine spline

      end module radial










