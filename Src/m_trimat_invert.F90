!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! General routine for inverting a matrix of arbitrary sizes.
! The algoritm follows: 10.1088/1749-4699/5/1/014009
! which contains an apparent general algoritm (not developed by the authors)
! 
! The idea is that we calculate in the following way (the order is important):
!   1. calculate Xn/Cn+1 and Yn/Bn-1.
!   2. calculate all Mnn
!   3. calculate all Mn-1n Mn+1n

! The algorithm has been developed so that it does not need any-more
! memory.

! It has fully been developed by Nick Papior Andersen, 2013
! Please contact the author of the code: nickpapior@gmail.com
! before use elsewhere!!!

! The routine can easily be converted to use another precision.

! When called you can specify two optional arguments which tells
! the module which parts you wish to calculate the inverse of.
! I.e. 
!   call invert_TriMat(M,Minv,sPart=2,ePart=3) 
! where parts >= 3. This will leverage a lot of computations,
! especially if parts >> 1 and only part of the matrix is wanted.


! After the inverted matrix has been calculated the input matrix
! looses the values in the diagonal
! The rest of the information is still retained
! If one only requests a part of the matrix, we retain some of the values
!   1. the lower tri-parts of the inverted tri-matrix contains
!       Xn/Cn+1 from [ePart+1:Parts]
!   2. the upper tri-parts of the inverted tri-matrix contains
!       Yn/Bn-1 from [1:sPart-1]
!   3. the lower tri-parts of the tri-matrix to be inverted contains
!       Xn/Cn+1 from [sPart:ePart]
!   4. the upper tri-parts of the tri-matrix to be inverted contains
!       Yn/Bn-1 from [sPart:ePart]
!   5. the A(1:sPart-1,1:sPart-1)         is retained IFF sPart > 1 
!   6. the A(ePart+1:Parts,ePart+1:Parts) is retained IFF ePart < Parts

module m_trimat_invert
  
  use class_zTriMat
  use precision, only: dp
  use m_pivot_array, only : Npiv, ipiv, init_pivot
  use m_pivot_array, only : clear_TriMat_inversion => clear_pivot 

  implicit none

  private
  private :: dp

  ! Used for BLAS calls (local variables)
  complex(dp), private, parameter :: z0 = cmplx( 0._dp, 0._dp, dp)
  complex(dp), private, parameter :: z1 = cmplx( 1._dp, 0._dp, dp)
  complex(dp), private, parameter :: zm1 = cmplx(-1._dp, 0._dp, dp)
#ifdef USE_GEMM3M
# define GEMM zgemm3m
#else
# define GEMM zgemm
#endif

  public :: invert_TriMat
  public :: init_TriMat_inversion
  public :: clear_TriMat_inversion

  ! For those inclined to do other things
  ! than simply inverting the matrix... :)
  public :: calc_Xn_div_Cn_p1
  public :: calc_Yn_div_Bn_m1
  public :: Xn_div_Cn_p1
  public :: Yn_div_Bn_m1
  public :: calc_Mnn_inv
  public :: calc_Mnm1n_inv
  public :: calc_Mnp1n_inv

contains

  subroutine invert_TriMat(M,Minv,calc_parts)
    type(zTriMat), intent(inout) :: M, Minv
    logical, intent(in), optional :: calc_parts(:)
    complex(dp), pointer :: Mpinv(:)
    integer :: lsPart, lePart
    integer :: sNm1, sNp1, n
    logical :: piv_initialized
    logical, allocatable :: lc_parts(:)

#ifndef TS_NOCHECKS
    if ( parts(M) /= parts(Minv) ) then
      call die('Could not calculate the inverse on non equal sized &
          &matrices')
    end if
    if ( parts(M) == 1 ) then
      call die('This matrix is not tri-diagonal')
    end if
    piv_initialized = .true.
    do n = 1 , parts(M)
      if ( Npiv < nrows_g(M,n) ) piv_initialized = .false.
    end do
    if ( .not. piv_initialized ) then
      call die('Pivoting array for inverting matrix not set.')
    end if
#endif

    ! Figure out if the calc_parts is correctly sized
    allocate(lc_parts(parts(M)))
    if ( present(calc_parts) ) then
#ifndef TS_NOCHECKS
      if ( size(calc_parts) /= parts(M) ) then
        call die('Wrong size of calculation parts. Please correct code')
      end if
#endif
      ! Copy over values
      lc_parts(:) = calc_parts(:)
      do n = 1 , parts(M)
        if ( lc_parts(n) ) then
          lsPart = n
          exit
        end if
      end do
      do n = parts(M) , 1 , -1
        if ( lc_parts(n) ) then
          lePart = n
          exit
        end if
      end do
    else
      lc_parts(:) = .true.
      lsPart = 1
      lePart = parts(M)
    end if

    call timer('TM_inv',1)

    ! Calculate all Xn/Cn+1
    do n = parts(M) - 1 , lsPart , -1 
      Mpinv => val(Minv,n+1,n+1)
      sNp1 = nrows_g(M,n+1)
      call calc_Xn_div_Cn_p1(M,Minv, n, Mpinv, sNp1**2 )
    end do
    ! Calculate all Yn/Bn-1
    do n = 2 , lePart
      Mpinv => val(Minv,n-1,n-1)
      sNm1 = nrows_g(M,n-1)
      call calc_Yn_div_Bn_m1(M,Minv, n, Mpinv, sNm1**2 )
    end do
    
    ! We calculate all Mnn
    ! Here it is permissable to overwrite the old A

    do n = lsPart , lePart
      if ( lc_parts(n) ) then
        call calc_Mnn_inv(M,Minv,n)
      end if
    end do

    ! ************ We have now calculated all diagonal parts of the 
    ! tri-diagonal matrix... **************************************

    do n = lsPart + 1 , lePart
      if ( lc_parts(n) ) then
        call calc_Mnm1n_inv(M,Minv,n)
      end if
    end do
    do n = lePart - 1 , lsPart , -1
      if ( lc_parts(n) ) then
        call calc_Mnp1n_inv(M,Minv,n)
      end if
    end do

    ! De-allocate variable to track calculated parts
    deallocate(lc_parts)

    call timer('TM_inv',2)
       
  end subroutine invert_TriMat

  ! We will calculate the diagonal inverse of the matrix.
  ! Before calling this routine the following routines *have*
  ! to be called:
  !    calc_Yn_div_Bn_m1(..., n=n, ...)
  !    calc_Xn_div_Cn_p1(..., n=n, ...)
  ! Input matrix M
  !    diagonal is lost, M[n,n] is *random*
  ! Output matrix Minv
  !    Minv[n,n] = M^(-1)[n,n]
  subroutine calc_Mnn_inv(M,Minv,n)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n
    ! Local variables
    complex(dp), pointer :: Mp(:), Mpinv(:)
    complex(dp), pointer :: Xn(:), Yn(:), Cn(:), Bn(:)
    integer :: sNm1, sN, sNp1, i

    if ( 1 < n ) sNm1 = nrows_g(M,n-1)
    sN = nrows_g(M,n)
    if ( n < parts(M) ) sNp1 = nrows_g(M,n+1)
    
    ! Retrieve Ann
    Mp => val(M,n,n)
    if ( n == 1 ) then
      ! First we calculate M11^-1
      ! Retrieve the X1/C2 array
      Xn => Xn_div_Cn_p1(Minv,n)
      ! The C2 array
      Cn => val(M,n,n+1)
      ! Calculate: A1 - X1
      call GEMM ('N','N',sN,sN,sNp1,zm1, &
          Cn,sN,Xn,sNp1,z1,Mp,sN)

    else if ( n == parts(M) ) then
      ! Retrieve the Yn/Bn-1 array
      Yn => Yn_div_Bn_m1(Minv,n)
      ! The Bn-1 array
      Bn => val(M,n,n-1)
      ! Calculate: An - Yn
      call GEMM ('N','N',sN,sN,sNm1,zm1, &
          Bn,sN,Yn,sNm1,z1,Mp,sN)

    else
      ! Retrieve the Xn/Cn+1 array
      Xn => Xn_div_Cn_p1(Minv,n)
      ! The Cn+1 array
      Cn => val(M,n,n+1)
      ! Calculate: An - Xn
      call GEMM ('N','N',sN,sN,sNp1,zm1, &
          Cn,sN,Xn,sNp1,z1,Mp,sN)
      ! Retrieve the Yn/Bn-1 array
      Yn => Yn_div_Bn_m1(Minv,n)
      ! The Bn-1 array
      Bn => val(M,n,n-1)
      ! Calculate: An - Xn - Yn
      call GEMM ('N','N',sN,sN,sNm1,zm1, &
          Bn,sN,Yn,sNm1,z1,Mp,sN)
  
    end if

    ! Retrive the position in the inverted matrix
    Mpinv => val(Minv,n,n)
    call zcopy(sN*sN, Mp, 1, Mpinv, 1)
    call zgetrf(sN, sN, Mpinv, sN, ipiv, i)
    if ( i == 0 ) then
      call zgetri(sN, Mpinv, sN, ipiv, Mp, sN*sN, i)
    else
      call die('Error on LU factorization of Mnn')
    end if
    if ( i /= 0 ) call die('Error on inverting Mnn')

  end subroutine calc_Mnn_inv

  ! Calculate the off-diagonal inverse of the matrix.
  ! Before calling this routine the following routines *have*
  ! to be called:
  !    calc_Mnn_inv(..., n=n, ...)
  !    calc_Xn_div_Cn_p1(..., n=n, ...)
  ! Input matrix M
  !    M[n+1,n] is destroyed
  ! Output matrix Minv
  !    Minv[n+1,n] = M^(-1)[n-1,n]
  subroutine calc_Mnp1n_inv(M,Minv,n)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n
    ! Local variables
    complex(dp), pointer :: Mp(:), Mpinv(:)
    complex(dp), pointer :: Xn(:)
    integer :: sN, sNp1

    if ( n < parts(M) ) then
      sNp1 = nrows_g(M,n+1)
    else
      ! We can/shall not calculate this
      return
    end if
    sN = nrows_g(M,n)

    ! *** we will now calculate Mn+1,n
    ! Copy over Xn/Cn+1
    Xn => Xn_div_Cn_p1(M,n)
    Mpinv => Xn_div_Cn_p1(Minv,n)

    call zcopy(sN*sNp1,Mpinv,1,Xn,1)

    ! Do matrix-multiplication
    Mp => val(Minv,n,n)
    ! Calculate: Xn/Cn+1 * Mnn
    call GEMM ('N','N',sNp1,sN,sN,zm1, &
        Xn,sNp1,Mp,sN,z0,Mpinv,sNp1)

  end subroutine calc_Mnp1n_inv

  ! Calculate the off-diagonal inverse of the matrix.
  ! Before calling this routine the following routines *have*
  ! to be called:
  !    calc_Mnn_inv(..., n=n, ...)
  !    calc_Yn_div_Bn_m1(..., n=n, ...)
  ! Input matrix M
  !    M[n-1,n] is destroyed
  ! Output matrix Minv
  !    Minv[n-1,n] = M^(-1)[n-1,n]
  subroutine calc_Mnm1n_inv(M,Minv,n)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n
    ! Local variables
    complex(dp), pointer :: Mp(:), Mpinv(:)
    complex(dp), pointer :: Yn(:)
    integer :: sN, sNm1

    if ( 1 < n ) then
      sNm1 = nrows_g(M,n-1)
    else
      ! We can/shall not calculate this
      return
    end if
    sN = nrows_g(M,n)

    ! Copy over Yn/Bn-1
    Yn => Yn_div_Bn_m1(M,n)
    Mpinv => Yn_div_Bn_m1(Minv,n)

    call zcopy(sN*sNm1,Mpinv,1,Yn,1)

    ! Do matrix-multiplication
    Mp => val(Minv,n,n)
    ! Calculate: Yn/Bn-1 * Mnn
    call GEMM ('N','N',sNm1,sN,sN,zm1, &
        Yn,sNm1,Mp,sN,z0,Mpinv,sNm1)
          
  end subroutine calc_Mnm1n_inv
    


  ! We will calculate the Xn/Cn+1 component of the 
  ! tri-diagonal inversion algorithm.
  ! The Xn/Cn+1 will be saved in the Minv n,n-1 (as that has
  ! the same size).
  ! Input matrix M
  !    nothing is changed
  ! Output matrix Minv
  !    Minv[n+1,n] = Xn/Cn+1
  ! zwork is altered (zwork may be Minv[n+1,n+1])
  subroutine calc_Xn_div_Cn_p1(M,Minv,n,zwork,nz)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, nz
    complex(dp), intent(inout) :: zwork(nz)
    ! Local variables
    complex(dp), pointer :: ztmp(:), Xn(:), Cnp2(:)
    integer :: sN, sNp1, sNp1SQ, sNp2, ierr
    character(len=50) :: cerr

#ifndef TS_NOCHECKS
    if ( n < 1 .or. parts(M) <= n .or. parts(M) /= parts(Minv) ) then
      call die('Could not calculate Xn on these matrices')
    end if
#endif
    ! Collect all matrix sizes for this step...
    sN = nrows_g(M,n)
    sNp1 = nrows_g(M,n+1)
    sNp1SQ = sNp1 ** 2
#ifndef TS_NOCHECKS
    if ( nz < sNp1SQ ) then
      call die('Work array in Xn calculation not sufficiently &
          &big.')
    end if
#endif

    ! Copy over the An+1 array
    ztmp => val(M, n+1, n+1)
    call zcopy(sNp1SQ, ztmp(1), 1, zwork(1), 1)

    ! If we should calculate X_N-1 then X_N == 0
    if ( n < parts(M) - 1 ) then
       ! Size...
       sNp2 =  nrows_g(M,n+2)
       ! Retrieve the Xn+1/Cn+2 array
       ztmp => Xn_div_Cn_p1(Minv,n+1)
       ! Retrieve the Cn+2 array
       Cnp2 => val(M,n+1,n+2)
       ! Calculate: An+1 - Xn+1
       call GEMM ('N','N',sNp1,sNp1,sNp2,zm1, &
           Cnp2,sNp1,ztmp,sNp2,z1,zwork,sNp1)
    end if

    ! Calculate Xn/Cn+1
    Xn => Xn_div_Cn_p1(Minv, n)
    ! Retrieve Cn+1 (stored in Cnp2 variable)
    Cnp2 => val(M, n+1, n)

    ! Split calculation into *fast*-slow
    ! For n,m n<m*2 zgesv is faster than zgetrf+zgetri+zgemm
    ! For n,m n>m*2 zgetrf+zgetri+zgemm is faster than zgesv
    if ( sNp1 * 2 > sN ) then
      ! Copy over the Bn array
      call zcopy(sN*sNp1, Cnp2(1), 1, Xn(1), 1)
      call zgesv(sNp1,sN,zwork,sNp1,ipiv,Xn,sNp1,ierr)
    else
      ! Here we know that sNp1 * 2 < sN and thus we
      ! can use the result array (Xn) as work array
      ! Since Xn has dimensions sNp1,sN > sNp1,sNp1
      call zgetrf(sNp1, sNp1, zwork, sNp1, ipiv, ierr)
      if ( ierr == 0 ) then
        call zgetri(sNp1, zwork, sNp1, ipiv, Xn, sNp1SQ, ierr)
        call GEMM ('N', 'N', sNp1, sN, sNp1, z1, &
            zwork, sNp1, Cnp2, sNp1, z0, Xn, sNp1)
      end if
    end if
    if ( ierr /= 0 ) then
      write(cerr,'(3(a,i0))') &
          'Error on inverting X',n,'/C',n+1,' with error: ',ierr
      call die(trim(cerr))
    end if

  end subroutine calc_Xn_div_Cn_p1

  function Xn_div_Cn_p1(M,n) result(Xn)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n
    complex(dp), pointer :: Xn(:)
    Xn => val(M,n+1,n)
  end function Xn_div_Cn_p1

  ! We will calculate the Yn/Bn-1 component of the 
  ! tri-diagonal inversion algorithm.
  ! The Yn/Bn-1 will be saved in the Minv n-1,n (as that has
  ! the same size).
  ! Input matrix M
  !    nothing is changed
  ! Output matrix Minv
  !    Minv[n-1,n] = Yn/Bn-1
  ! zwork is altered (zwork may be Minv[n-1,n-1])
  subroutine calc_Yn_div_Bn_m1(M,Minv,n,zwork,nz)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, nz
    complex(dp), intent(inout) :: zwork(nz)
    ! Local variables
    complex(dp), pointer :: ztmp(:), Yn(:), Bnm2(:)
    integer :: sN, sNm1, sNm1SQ, sNm2, ierr
    character(len=50) :: cerr

#ifndef TS_NOCHECKS
    if ( n < 2 .or. parts(M) < n .or. parts(M) /= parts(Minv) ) then
      call die('Could not calculate Yn on these matrices')
    end if
#endif
    ! Collect all matrix sizes for this step...
    sN = nrows_g(M,n)
    sNm1 = nrows_g(M,n-1)
    sNm1SQ = sNm1 ** 2
#ifndef TS_NOCHECKS
    if ( nz < sNm1SQ ) then
      call die('Work array in Yn calculation not sufficiently &
          &big.')
    end if
#endif

    ! Copy over the An-1 array
    ztmp => val(M, n-1, n-1)
    call zcopy(sNm1SQ, ztmp(1), 1, zwork(1), 1)

    if ( 2 < n ) then
      ! Size...
      sNm2 =  nrows_g(M,n-2)
      ! Retrieve the Yn-1/Bn-2 array
      ztmp => Yn_div_Bn_m1(Minv,n-1)
      ! Retrieve the Bn-2 array
      Bnm2 => val(M,n-1,n-2)
      ! Calculate: An-1 - Yn-1
      call GEMM ('N','N',sNm1,sNm1,sNm2,zm1, &
          Bnm2,sNm1,ztmp,sNm2,z1,zwork,sNm1)
    end if

    ! Calculate Yn/Bn-1
    Yn => Yn_div_Bn_m1(Minv, n)
    ! Retrieve Bn-1 (stored in Bnm2 variable)
    Bnm2 => val(M, n-1, n)

    ! Split calculation into *fast*-slow
    ! For n,m n<m*2 zgesv is faster than zgetrf+zgetri+zgemm
    ! For n,m n>m*2 zgetrf+zgetri+zgemm is faster than zgesv
    if ( sNm1 * 2 > sN ) then
      call zcopy(sN*sNm1, Bnm2(1), 1, Yn(1), 1)
      call zgesv(sNm1,sN,zwork,sNm1,ipiv,Yn,sNm1,ierr)
    else
      ! Here we know that sNm1 * 2 < sN and thus we
      ! can use the result array (Yn) as work array
      ! Since Yn has dimensions sNm1,sN > sNm1,sNm1
      call zgetrf(sNm1, sNm1, zwork, sNm1, ipiv, ierr)
      if ( ierr == 0 ) then
        call zgetri(sNm1, zwork, sNm1, ipiv, Yn, sNm1SQ, ierr)
        call GEMM ('N', 'N', sNm1, sN, sNm1, z1, &
            zwork, sNm1, Bnm2, sNm1, z0, Yn, sNm1)
      end if
    end if
    if ( ierr /= 0 ) then
      write(cerr,'(3(a,i0))') &
          'Error on inverting Y',n,'/B',n-1,' with error: ',ierr
      call die(trim(cerr))
    end if

  end subroutine calc_Yn_div_Bn_m1

  function Yn_div_Bn_m1(M,n) result(Yn)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n
    complex(dp), pointer :: Yn(:)
    Yn => val(M,n-1,n)
  end function Yn_div_Bn_m1

  ! We initialize the pivoting array for rotating the inversion
  subroutine init_TriMat_inversion(M)
    type(zTriMat), intent(in) :: M
    integer :: i, N

    N = nrows_g(M,1)
    do i = 2 , parts(M)
      N = max(N, nrows_g(M,i))
    end do

    call init_pivot(N)

  end subroutine init_TriMat_inversion

#undef GEMM

end module m_trimat_invert
    
