subroutine filter( l, nr, r, f, kmax, norm_opt,emin )

! Filters out high-k components of a radial function
! J.M.Soler. July 2005
!------------------------------------------------------------
! Input:
!  integer l     : Angular momentum of the function to filter
!  integer nr    : Number of radial points
!  real*8  r(nr) : Radial mesh, in ascending order
!  real*8  kmax  : Cutoff in reciprocal space
!  integer norm_opt : Option to renormalize after filtering:
!                       0 => do not renormalize
!                       1 => renormalize f
!                       2 => renormalize f**2
!  real*8  emin_in: Min. norm. of filt. eigv.
! Input/output:
!  real*8  f(nr) : Radial function to be filtered
!------------------------------------------------------------

!  use plot_module

  implicit none
  integer, parameter :: dp = kind(1.d0)

! Arguments
  integer,  intent(in)    :: l
  integer,  intent(in)    :: nr
  real(dp), intent(in)    :: kmax
  real(dp), intent(in)    :: r(nr)
  real(dp), intent(inout) :: f(nr)
  integer,  intent(in)    :: norm_opt
  real(dp), intent(in)    :: emin

! Internal parameters
  real(dp), parameter :: pkr  = 1.00_dp ! Num. polyn. / (kmax*rmax)
  real(dp), parameter :: emin_default = 0.999_dp ! Min. norm. of filt. eigv.

! Internal variables and arrays
  integer :: i, j, ir, ix, jx, lp, m, n
  real(dp):: fnorm, f0norm, krmax, p0, pi, rmax, y,aa,aamax
  integer, allocatable:: indx(:)
  real(dp),allocatable:: a(:,:), ar(:,:), aux(:), e(:), &
                         f0(:), p(:,:), wx(:), x(:)

! Debugging variables
!  real(dp):: aa, aamax, af, axy, da, damax, df, &
!             fa, ff0(nr), fmax, fmin , fnorm

! Interfaces to external functions
!   interface
!     DOUBLE PRECISION FUNCTION BESSPH (L,X)
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!       integer  l
!     END
!     DOUBLE PRECISION FUNCTION PLGNDR(L,M,X)
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     END
!     SUBROUTINE XLGNDR (X,W,N)
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!       DOUBLE PRECISION X(N),W(N)
!     END
!   end interface
  real(dp),external :: bessph,plgndr
  
! Fix order of filter kernel by empirical rule
  rmax = r(nr)
  krmax = kmax * rmax
  n = l/2 + 1 + nint(krmax*pkr)
  !write(6,*) "n,krmax,rmax:",n,krmax,rmax

! Allocate arrays
  allocate( a(n,n), ar(nr,n), aux(n), e(n), f0(n), indx(n), &
            p(n,n), wx(n), x(n) )

! Find special points and weights for Gauss-Legendre integration
  lp = mod(l,2)
  call xlgndr( x, wx, 2*n-lp )

! Calculate the "filter" kernel at the special points
  do ix = 1,n
    do jx = 1,ix
      y = krmax * x(ix) * x(jx)
      a(ix,jx) = y * bessph( l, y )
      a(jx,ix) = a(ix,jx)
    end do
  end do

! Find normalized Legendre polynomials at special points (times weight)
  do i = 1,n
    j = 2*i - lp - 1
    p0 = sqrt( real( 2*j+1, dp ) )
    do ix = 1,n
      p(ix,i) = p0 * plgndr( j, 0, x(ix) ) * wx(ix)
    end do
  end do

! Expand the kernel in orthonormal Legendre polynomials: pT*a*p
  a = matmul( transpose(p), matmul(a,p) )

! Remove integration weight from Legendre polynomials
  forall (i=1:n, ix=1:n) p(ix,i) = p(ix,i) / wx(ix)

! Diagonalize the filter kernel
  call tred2( a, n, n, e, aux )
  call tqli( e, aux, n, n, a )

! Order eigenvalues and eigenvectors with decreasing eigval**2
  call ordix( -e**2, 1, n, indx )
  call order( e, 1, n, indx )
  call order( a, n, n, indx )

! Print eigenvalues for debugging
  pi = acos(-1._dp)
  !print'(a,i4,2f18.12)', &
  !  ('filter: i,e,e2=',i,sqrt(2*krmax/pi)*e(i),2*krmax/pi*e(i)**2,i=1,n)

! Renormalize eigenvalues so that e=1 => perfect confinement 
  pi = acos(-1._dp)
  e = 2*krmax/pi * e**2

! Find how many eigenvalues are within filter tolerance
  m = 0
  do i = 1,n
    if (e(i)<emin) exit
    m = i
  end do
!  print*, 'filter: n, m =', n, m

! Find Legendre polynomials at mesh points
  ar = 0.0_dp
  do i = 1,n
    j = 2*i - lp - 1
    p0 = sqrt( real( 2*j+1, dp ) )
    do ir = 1,nr
      ar(ir,i) = ar(ir,i) + p0 * plgndr( j, 0, r(ir)/rmax )
    end do
  end do

! Find eigenvectors at mesh points
  ar = matmul( ar, a )

! Find eigenvectors at special points
  a = matmul( p, a )

! Check that eigenvectors are orthonormal for debugging
 aamax = 0
  do i = 1,n
    do j = 1,i
!!      aa = sum( wx(:) * p(:,i) * p(:,j) )
      aa = sum( wx(:) * a(:,i) * a(:,j) )
!!      print'(a,2i6,f15.9)', 'filter: i,j,<i|j>=', i, j, aa
      if (i==j) aa = aa-1
      aamax = max( aamax, abs(aa) )
    end do
  end do
  !print*, 'filter: n, aamax =', n, aamax

! Plot first eigenvectors
!  call window( 0._dp, 1._dp, -3._dp, 3._dp, 0 )
!  call axes( 0._dp, 0.1_dp, 0._dp, 0.5_dp )
!  do i = 1,min(n,5)
!!    call plot( n, x, a(:,i) )
!    call plot( nr, r/rmax, ar(:,i)*sign(1._dp,ar(nr,i)) )
!!    call plot( nr-1, r(2:nr)/rmax, ar(2:nr,i)/(r(2:nr)/rmax) )
!  end do
!  call show(0)

! Find function to be filtered at special points
  f0(:) = interpolation( f, r/rmax, x )

!  print'(a,2f12.6)', ('filter: x,f0=',x(i),f0(i),i=1,n)

! Filter radial function (times x)
!  ff0 = f
  f = 0
  do i = 1,m
    f(1:nr) = f(1:nr) + ar(1:nr,i) * sum(wx*x*f0*a(:,i))
  end do

! Renormalize filtered function
  if (norm_opt==1) then
    f0norm = sum( matmul(wx*x*f0,a(:,1:n)) )
    fnorm  = sum( matmul(wx*x*f0,a(:,1:m)) )
    f = f * f0norm/fnorm
  else if (norm_opt==2) then
    f0norm = sum( matmul(wx*x*f0,a(:,1:n))**2 )
    fnorm  = sum( matmul(wx*x*f0,a(:,1:m))**2 )
    f = f * sqrt(f0norm/fnorm)
  else if (norm_opt/=0) then
    stop 'filter: ERROR: invalid value of norm_opt'
  end if

!  if (norm_opt==1 .or. norm_opt==2) &
!    print*, 'filter: f0norm, fnorm =', f0norm*rmax**3, fnorm*rmax**3

! Divide filtered function (times x) by x=r/rmax
  f(2:nr) = f(2:nr) / (r(2:nr)/rmax)
  if (r(1)==0._dp) then
    if (l==0) then
      ! Linear extrapolation to r=0
      f(1) = (r(3)*f(2) - r(2)*f(3)) / (r(3) - r(2))
    else
      f(1) = 0
    end if
  else
    f(1) = f(1) / (r(1)/rmax)
  end if



!  print*, 'filter: dfmean =', sqrt(sum((f-ff0)**2)/nr)

! Plot original and filtered functions for debugging
!  fmin = minval(ff0)
!  fmax = maxval(ff0)
!  df = fmax - fmin

!  call window( 0._dp, 1._dp, fmin-df/10, fmax+df/10, 0 )
!  call axes( 0._dp, 0.1_dp, 0._dp, df/10 )
!!  call plot( n, x, f0 )
!  call plot( nr, r/rmax, ff0, f )
!!  call plot( nr-1, r(2:nr)/rmax, ff0(2:nr), f(2:nr) )
!  call show(0)

! Deallocate arrays
  deallocate( a, ar, aux, e, f0, indx, p, wx, x )

contains

  function interpolation( yold, xold, xnew ) result(ynew)

    ! Interpolates function yold(xold) at new points xnew
    ! Assumes that xold and xnew are in increasing order

    implicit none
    real(dp), intent(in) :: yold(:), xold(:), xnew(:)
    real(dp) :: ynew(size(xnew))

    integer :: i1, i2, inew, iold, nnew, nold

    nold = size(xold)
    nnew = size(xnew)

    iold = 1
    do inew = 1,nnew
      ! Find iold such that xold(iold) < xnew(inew) < xold(iold+1)
      do
        if (iold+1==nold .or. xold(iold+1) > xnew(inew)) exit
        iold = iold + 1
      end do
      ! Four-point Lagrange interpolation (3-point at extreme intervals)
      i1 = max(iold-1,1)
      i2 = min(iold+2,nold)
      ynew(inew) = lagrange( yold(i1:i2), xold(i1:i2), xnew(inew) )
    end do

  end function interpolation

  function lagrange( yold, xold, xnew ) result(ynew)

    ! Lagrange interpolation of function yold(xold) at point xnew
    ! Based on routine polint of Numerical Recipes

    implicit none
    real(dp), intent(in) :: yold(:), xold(:), xnew
    real(dp) :: ynew

    integer :: i, i0, m, n
    real(dp):: c(size(xold)), d(size(xold)), den, dy, ho, hp, w

    n = size(xold)
    c = yold
    d = yold
    i0 = n/2
    ynew = yold(i0)
    i0 = i0 - 1
    do m = 1,n-1
      do i=1,n-m
        ho = xold(i) - xnew
        hp = xold(i+m) - xnew
        w = c(i+1) - d(i)
        den = ho - hp
        if (den==0._dp) stop 'filter:lagrange: ERROR: den=0'
        d(i) = hp * w / den
        c(i) = ho * w / den
      end do
      if (2*i0 < n-m) then
        dy = c(i0+1)
      else
        dy = d(i0)
        i0 = i0 - 1
      endif
      ynew = ynew + dy
    end do

  end function lagrange

end subroutine filter
