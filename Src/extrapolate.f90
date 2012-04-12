!!@LICENSE
!
! *******************************************************************
! SUBROUTINE extrapolate( na, n, cell, xa, cell0, x0, c )
! *******************************************************************
! Finds optimal coefficients for an approximate expasion of the form
! D_0 = sum_i c_i D_i, where D_i is the density matrix in the i'th
! previous iteration, or any other function of the atomic positions. 
! In practice, given points x_0 and x_i, i=1,...,n, it finds the
! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
! sum_i c_i = 1. Unit cell vectors are included for completeness of
! the geometry specification only. Though not used in this version,
! this routine can be used if cell vectors change (see algorithms).
! Written by J.M.Soler. Feb.2010.
! ************************* INPUT ***********************************
!  integer  na          ! Number of atoms
!  integer  n           ! Number of previous atomic positions
!  real(dp) cell(3,3,n) ! n previous unit cell vectors
!  real(dp) xa(3,na,n)  ! n previous atomic positions (most recent first)
!  real(dp) cell0(3,3)  ! Present unit cell vectors
!  real(dp) x0(3,na)    ! Present atomic positions
! ************************* OUTPUT **********************************
!  real(dp) c(n)        ! Expansion coefficients
! ************************ UNITS ************************************
! Unit of distance is arbitrary, but must be the same in all arguments
! ********* BEHAVIOUR ***********************************************
! - Returns without any action if n<1 or na<1
! - Stops with an error message if matrix inversion fails
! ********* DEPENDENCIES ********************************************
! Routines called: 
!   inver     : Matrix inversion
! Modules used:
!   precision : defines parameter 'dp' (double precision real kind)
!   sys       : provides the stopping subroutine 'die'
! ********* ALGORITHMS **********************************************
! Using a linear approximation for D(x), imposing that sum(c)=1, and
! assuming that <dD/dx_i*dD/dx_j>=0, <(dD/dx_i)**2>=<(dD/dx_j)**2>
! (i.e. assuming no knowledge on the parcial derivatives), we find
! (D(x0)-D(sum_i c_i*x_i))**2 = const * (x0-sum_i c_i*x_i)**2
! Therefore, given points x_0 and x_i, i=1,...,n, we look for the
! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
! sum_i c_i = 1. It is straighforward to show that this reduces to 
! solving S*c=s0, where S_ij=x_i*x_j and s0_i=x_i*x0, under the
! constraint sum(c)=1. To impose it, we rewrite the expansion as
! xmean + sum_i c_i*(x_i-xmean), where xmean=(1/n)*sum_i x_i.
! Since the vectors (x_i-xmean) are linearly dependent, the matrix
! S_ij=(x_i-xmean)*(x_j-xmean) is singular. Therefore, we substitute
! the last row of S and s0 by 1, thus imposing sum(c)=1.
! Unit cell vectors are not used in this version, but this routine
! can be safely used even if cell vectors change, since D depends
! on the interatomic distances, and it can be shown that
! sum_ia,ja (xa(:,ia,i)-xa(:,ja,i))*(xa(:,ia,j)-xa(:,ja,j))
!   = 2*na*sum_ia (xa(:,ia,i)-xmean(:,ia))*(xa(:,ia,j)-xmean(:,ia))
!   = 2*na*S_ij
! This implies that approximating the present ineratomic distances,
! as an expansion of previous ones, is equivalent to approximating
! the atomic positions.
! *******************************************************************

SUBROUTINE extrapolate( na, n, cell, xa, cell0, x0, c )

! Used procedures and parameters
  USE sys,       only: die           ! Termination routine
  USE precision, only: dp            ! Double precision real kind

! Passed arguments
  implicit none
  integer, intent(in) :: na          ! Number of atoms
  integer, intent(in) :: n           ! Number of previous atomic positions
  real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
  real(dp),intent(in) :: xa(3,na,n)  ! n previous atomic positions
  real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
  real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
  real(dp),intent(out):: c(n)        ! Expansion coefficients

! Internal variables and arrays
  integer :: i, ierr, j, m
  real(dp):: s(n,n), s0(n), si(n,n), xmean(3,na)

! Trap special cases
  if (na<1 .or. n<1) then
    return
  else if (n==1) then
    c(1) = 1
    return
  end if

! Find average of previous positions
  xmean = sum(xa,dim=3) / n

! Find matrix s of dot products. Subtract xmean to place origin within the
! hyperplane of x vectors
  do j = 1,n
    s0(j) = sum( (x0-xmean) * (xa(:,:,j)-xmean) )
    do i = j,n
      s(i,j) = sum( (xa(:,:,i)-xmean) * (xa(:,:,j)-xmean) )
      s(j,i) = s(i,j)
    end do
  end do

! Find the largest number of (first) m linearly independent vectors xa-xmean.
! Notice that m<n because we have subtracted xmean. Optimally, we should not
! restrict ourselves to the first (most recent) m vectors, but that would
! complicate the code too much.
  do m = n-1,0,-1
    if (m==0) exit                   ! Do not try to invert a 0x0 'matrix'
    call inver( s, si, m, n, ierr )
    if (ierr==0) exit                ! Not too elegant, but it will do.
  end do

! Trap the case in which the first two xa vectors are equal
  if (m==0) then  ! Just use most recent point only
    c(1) = 1
    c(2:n) = 0
    return
  end if

! Set one more row for equation sum(c)=1.
  m = m+1
  s0(m) = 1
  s(m,1:m) = 1

! Invert the full equations matrix
  call inver( s, si, m, n, ierr )
  if (ierr/=0) call die('extrapolate: ERROR: matrix inversion failed')

! Find expansion coefficients
  c(1:m) = matmul( si(1:m,1:m), s0(1:m) )
  c(m+1:n) = 0

END SUBROUTINE extrapolate

