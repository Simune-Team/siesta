!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
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
! before use

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
!   3. the A(1,1)         is retained IFF sPart > 1 
!   4. the A(Parts,Parts) is retained IFF ePart < Parts

module m_trimat_invert
  
  use class_zTriMat
  use precision, only: dp

  implicit none

  private
  private :: dp

  ! Current size of the pivoting arrays
  integer, private, save          :: Npiv = 0
  ! The pivoting array
  integer, private, save, pointer :: ipiv(:) => null()

  ! Used for BLAS calls (local variables)
  complex(dp), private, parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), private, parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), private, parameter :: zm1 = dcmplx(-1._dp, 0._dp)

  public :: invert_TriMat
  public :: init_TriMat_inversion
  public :: clear_TriMat_inversion

contains

  subroutine invert_TriMat(M,Minv,sPart,ePart)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in), optional :: sPart, ePart
    complex(dp), pointer :: Mpinv(:), Mp(:), XYn(:), CB(:)
    integer :: lsPart, lePart
    integer :: sNm1, sNp1, n
    logical :: piv_initialized

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

    call timer('TM_inv',1)

    lsPart = 1
    if ( present(sPart) ) lsPart = sPart
    lePart = parts(M)
    if ( present(ePart) ) lePart = ePart

    ! Calculate all Xn/Cn+1
    do n = parts(M) - 1 , lsPart , -1 
       Mpinv => val(Minv,n+1,n+1)
       sNp1 = nrows_g(M,n+1)
       call calc_Xn_div_Cn_p1(M,Minv, n, Mpinv,sNp1**2 )
    end do
    ! Calculate all Yn/Bn-1
    do n = 2 , lePart
       Mpinv => val(Minv,n-1,n-1)
       sNm1 = nrows_g(M,n-1)
       call calc_Yn_div_Bn_m1(M,Minv, n, Mpinv, sNm1**2)
    end do
    
    ! We calculate all Mnn
    ! Here it is permissable to overwrite the old A

    do n = lsPart , lePart
       call calc_Mnn_inv(M,Minv,n)
    end do

    ! ************ We have now calculated all diagonal parts of the 
    ! tri-diagonal matrix... **************************************

    do n = lsPart + 1 , lePart
       call calc_Mnm1n_inv(M,Minv,n)
    end do
    do n = lsPart , lePart - 1
       call calc_Mnp1n_inv(M,Minv,n)
    end do

    call timer('TM_inv',2)
       
  end subroutine invert_TriMat

  subroutine calc_Mnn_inv(M,Minv,n)
    use m_mat_invert
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n
    ! Local variables
    complex(dp), pointer :: Mp(:), Mpinv(:)
    complex(dp), pointer :: Xn(:), Yn(:), Cn(:), Bn(:)
    integer :: sNm1, sN, sNp1, ierr

    if ( 1 < n )        sNm1 = nrows_g(M,n-1)
                        sN   = nrows_g(M,n)
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
       call zgemm('N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)

    else if ( n == parts(M) ) then

       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Yn
       call zgemm('N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)

    else
       ! Retrieve the Xn/Cn+1 array
       Xn => Xn_div_Cn_p1(Minv,n)
       ! The Cn+1 array
       Cn => val(M,n,n+1)
       ! Calculate: An - Xn
       call zgemm('N','N',sN,sN,sNp1, &
            zm1, Cn,sN, Xn,sNp1,z1, Mp,sN)
       ! Retrieve the Yn/Bn-1 array
       Yn => Yn_div_Bn_m1(Minv,n)
       ! The Bn-1 array
       Bn => val(M,n,n-1)
       ! Calculate: An - Xn - Yn
       call zgemm('N','N',sN,sN,sNm1, &
            zm1, Bn,sN, Yn,sNm1,z1, Mp,sN)
  
    end if

    ! Retrive the position in the inverted matrix
    Mpinv => val(Minv,n,n)
    call mat_invert(Mp,Mpinv,sN,method = MI_WORK )

  end subroutine calc_Mnn_inv

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
    sN    = nrows_g(M,n)

    ! *** we will now calculate Mn+1,n
    ! Copy over Xn/Cn+1
    Xn    => Xn_div_Cn_p1(M   ,n)
    Mpinv => Xn_div_Cn_p1(Minv,n)
    Xn(:) =  Mpinv(:)
    ! Do matrix-multiplication
    Mp    => val(Minv,n,n)
    ! Calculate: Xn/Cn+1 * Mnn
    call zgemm('N','N',sNp1,sN,sN, &
         zm1, Xn,sNp1, Mp,sN,z0, Mpinv,sNp1)

  end subroutine calc_Mnp1n_inv

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
    sN   = nrows_g(M,n)

    ! Copy over Yn/Bn-1
    Yn    => Yn_div_Bn_m1(M   ,n)
    Mpinv => Yn_div_Bn_m1(Minv,n)
    Yn(:) =  Mpinv(:)
    ! Do matrix-multiplication
    Mp     => val(Minv,n,n)
    ! Calculate: Yn/Bn-1 * Mnn
    call zgemm('N','N',sNm1,sN,sN, &
         zm1, Yn,sNm1, Mp,sN,z0, Mpinv,sNm1)
          
  end subroutine calc_Mnm1n_inv
    


  ! We will calculate the Xn/Cn+1 component of the 
  ! tri-diagonal inversion algorithm.
  ! The Xn/Cn+1 will be saved in the Minv n,n-1 (as that has
  ! the same size).
  subroutine calc_Xn_div_Cn_p1(M,Minv,n,zwork,nz)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, nz
    complex(dp), intent(inout) :: zwork(nz)
    ! Local variables
    complex(dp), pointer :: ztmp(:), Xn(:), Cnp2(:)
    integer :: sN, sNp1, sNp1SQ, sNp2, ierr

    if ( n < 1 .or. parts(M) <= n .or. parts(M) /= parts(Minv) ) then
       call die('Could not calculate Xn on these matrices')
    end if
    ! Collect all matrix sizes for this step...
    sN     = nrows_g(M,n)
    sNp1   = nrows_g(M,n+1)
    sNp1SQ = sNp1 ** 2
    if ( nz < sNp1SQ ) then
       call die('Work array in Xn calculation not sufficiently &
            &big.')
    end if

    ! Copy over the Bn array
    ztmp  => val(M   ,n+1,n)
    ! This is where the inverted matrix will be located 
    Xn    => Xn_div_Cn_p1(Minv,n)
    Xn(:) =  ztmp(:)

    ! Copy over the An+1 array
    ztmp            => val(M,n+1,n+1)
    zwork(1:sNp1SQ) =  ztmp(:)

    ! If we should calculate X_N-1 then X_N == 0
    if ( n < parts(M) - 1 ) then
       ! Size...
       sNp2 =  nrows_g(M,n+2)
       ! Retrieve the Xn+1/Cn+2 array
       ztmp => Xn_div_Cn_p1(Minv,n+1)
       ! Retrieve the Cn+2 array
       Cnp2 => val(M,n+1,n+2)
       ! Calculate: An+1 - Xn+1
       call zgemm('N','N',sNp1,sNp1,sNp2, &
            zm1, Cnp2,sNp1, ztmp,sNp2,z1, zwork,sNp1)
    end if

    ! Calculate Xn/Cn+1
    call zgesv(sNp1,sN,zwork,sNp1,ipiv,Xn,sNp1,ierr)
    if ( ierr /= 0 ) call die('Error on inverting Xn/Cn+1')

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
  subroutine calc_Yn_div_Bn_m1(M,Minv,n,zwork,nz)
    type(zTriMat), intent(inout) :: M, Minv
    integer, intent(in) :: n, nz
    complex(dp), intent(inout) :: zwork(nz)
    ! Local variables
    complex(dp), pointer :: ztmp(:), Yn(:), Bnm2(:)
    integer :: sN, sNm1, sNm1SQ, sNm2, ierr

    if ( n < 2 .or. parts(M) < n .or. parts(M) /= parts(Minv) ) then
       call die('Could not calculate Xn on these matrices')
    end if
    ! Collect all matrix sizes for this step...
    sN     = nrows_g(M,n)
    sNm1   = nrows_g(M,n-1)
    sNm1SQ = sNm1 ** 2
    if ( nz < sNm1SQ ) then
       call die('Work array in Xn calculation not sufficiently &
            &big.')
    end if

    ! Copy over the Cn array
    ztmp  => val(M   ,n-1,n)
    ! This is where the inverted matrix will be located 
    Yn    => Yn_div_Bn_m1(Minv,n)
    Yn(:) =  ztmp(:)

    ! Copy over the An-1 array
    ztmp            => val(M,n-1,n-1)
    zwork(1:sNm1SQ) =  ztmp(:)

    ! If we should calculate Y_2 then Y_1 == 0
    if ( 2 < n ) then
       ! Size...
       sNm2 =  nrows_g(M,n-2)
       ! Retrieve the Yn-1/Bn-2 array
       ztmp => Yn_div_Bn_m1(Minv,n-1)
       ! Retrieve the Bn-2 array
       Bnm2 => val(M,n-1,n-2)
       ! Calculate: An-1 - Yn-1
       call zgemm('N','N',sNm1,sNm1,sNm2, &
            zm1, Bnm2,sNm1, ztmp,sNm2,z1, zwork,sNm1)
    end if

    ! Calculate Yn/Bn-1
    call zgesv(sNm1,sN,zwork,sNm1,ipiv,Yn,sNm1,ierr)
    if ( ierr /= 0 ) call die('Error on inverting Yn/Bn-1')

  end subroutine calc_Yn_div_Bn_m1

  function Yn_div_Bn_m1(M,n) result(Yn)
    type(zTriMat), intent(in) :: M
    integer, intent(in) :: n
    complex(dp), pointer :: Yn(:)
    Yn => val(M,n-1,n)
  end function Yn_div_Bn_m1

  ! We initialize the pivoting array for rotating the inversion
  subroutine init_TriMat_inversion(M)
    use alloc, only : re_alloc
    type(zTriMat), intent(in) :: M
    integer :: i
    Npiv = 0
    do i = 1 , parts(M)
       if ( nrows_g(M,i) > Npiv ) then
          Npiv = nrows_g(M,i)
       end if
    end do

    ! Allocate space for the pivoting array
    call re_alloc(ipiv,1, Npiv, &
         name="TriMat_piv",routine='TriMatInversion')

  end subroutine init_TriMat_inversion

  subroutine clear_TriMat_inversion()
    use alloc, only: de_alloc
    Npiv = 0
    ! Deallocate the pivoting array
    call de_alloc(ipiv, &
         name="TriMat_piv",routine='TriMatInversion')
  end subroutine clear_TriMat_inversion

end module m_trimat_invert
    
