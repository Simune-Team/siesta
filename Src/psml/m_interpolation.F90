module m_interpolation

use external_interfaces, only: die

implicit none

!
!  Implements simple polynomial interpolation
!
public :: interpolate
private

integer, parameter :: dp = selected_real_kind(10,100)

CONTAINS

subroutine interpolate(x,y,r,val,debug)
real(dp), intent(in) :: x(:), y(:)
real(dp), intent(in) :: r
real(dp), intent(out):: val
logical, intent(in), optional :: debug

integer, save :: i0
integer :: npts, nmin, nmax, nn
integer, parameter :: npoint = 2  ! basis for interpolation order
real(dp)  :: dy

npts = size(x)
if (size(y) /= npts) call die("x and y not conformable in interpolate")

! Find closest point in x array so that x(i0) < r
call hunt(x,npts,r,i0)
!
! Choose points to interpolate over.
! At the extremes of the range we get quadratic interpolation
! The maximum order is 2*npoint + 1 (cuartic polynomial)
! It is in principle possible to extrapolate, but with reduced
! order (and all the dangers pertaining to extrapolation!)
!
nmin=max(1,i0-npoint)
nmax=min(npts,i0+npoint)
nn=nmax-nmin+1
!
call polint(x(nmin:),y(nmin:),nn,r,val,dy)

if (present(debug)) then
   if (debug) then
      print "(a,/,2g20.10,3i3,g20.10)", &
          "r ,r-x(i0), i0, nmin, nmax, d: ", &
           r, r-x(i0), i0, nmin, nmax, dy
   endif
endif
end subroutine interpolate

      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL(dp), intent(in) ::  x,xx(n)

      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
    END SUBROUTINE hunt
!  (C) Copr. 1986-92 Numerical Recipes Software

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY) 
!*****************************************************************
! Polynomial interpolation. Modified and adapted to double 
! precision from same routine of Numerical Recipes.
! D. Sanchez-Portal, Oct. 1996
!*****************************************************************
! Input:
!   real*8  XA(N) : x values of the function y(x) to interpolate
!   real*8  YA(N) : y values of the function y(x) to interpolate
!   integer N     : Number of data points
!   real*8  X     : x value at which the interpolation is desired
! Output:
!   real*8  Y     : interpolated value of y(x) at X
!   real*8  DY    : accuracy estimate
!*****************************************************************

      IMPLICIT NONE
      INTEGER  :: N
      REAL(dp) :: XA(N),YA(N), X, Y, DY

      INTEGER  :: I, M, NS
      REAL(dp) :: C(N), D(N), DEN, DIF, DIFT, HO, HP, W
      REAL(dp), PARAMETER :: ZERO=0.D0

      NS=1
      DIF=ABS(X-XA(1))
      DO I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
      END DO ! I
      Y=YA(NS)
      NS=NS-1
      DO M=1,N-1
        DO I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF (DEN.EQ.ZERO) call die('polint: ERROR. Two XAs are equal')
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
        END DO ! I
        IF (2*NS.LT.N-M) THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
      END DO ! M

      END SUBROUTINE POLINT

end module m_interpolation

#ifdef __TEST__
program test_interp
use m_interpolation
integer, parameter :: dp = selected_real_kind(10,100)
real(dp) :: x(7) = (/ ((i-1)*0.15, i=1,7) /)
real(dp) :: y(7)
real(dp) :: xi, vi, v
y = sin(x)
do
   print *, "Enter value x in [0,0.9]"
   read *, xi
   call interpolate(x,y,xi,vi,debug=.true.)
   v = sin(xi)
   print "(a,2f15.8,g20.10)", "Interp, true, diff: ", vi, v, (vi-v)
enddo
end program test_interp

      subroutine die(str)
      character(len=*), intent(in), optional :: str
      if (present(str)) then
         write(*,*) str
      endif
      STOP
      end subroutine die

#endif
