module m_mat_invert
  
  use precision, only : dp
  
  implicit none

  ! The maximum dimensionality of the problem before we turn to a direct inversion algorithm
  integer, private :: N_MAX = 40
  integer, private, save :: Npiv
  integer, private, save, pointer :: ipiv(:) => null()

  ! Used for BLAS calls (local variables)
  complex(dp), private, parameter :: z0  = dcmplx( 0._dp, 0._dp)
  complex(dp), private, parameter :: z1  = dcmplx( 1._dp, 0._dp)
  complex(dp), private, parameter :: zm1 = dcmplx(-1._dp, 0._dp)

contains

  recursive subroutine mat_invert(M, zwork, no)
    use intrinsic_missing, only : EYE
    complex(dp), pointer :: M(:), zwork(:)
    integer, intent(in) :: no ! Size of problem

    complex(dp), pointer :: A1(:), C2(:), B1(:), A2(:)
    complex(dp), pointer :: X1(:), Y2(:), t1(:), t2(:)
    integer :: sA1, eA1, sC2, eC2
    integer :: sB1, eB1, sA2, eA2
    integer :: n1, n2, i1, i2, idx, i, j

    ! M is the matrix we wish to invert using the tri-diagonal
    ! routine.
    ! This will enable us to half no in order to "faster" 
    ! invert the matrix

    if ( no <= N_MAX ) then
       idx = 0
       do j = 1 , no 
          do i = 1 , no 
             idx = idx + 1
             zwork(idx) = M(idx)
             if ( i == j ) then
                M(idx) = z1
             else
                M(idx) = z0
             end if
          end do
       end do
       call zgesv(no,no,zwork,no,ipiv,M,no,i)
       if ( i /= 0 ) call die('Error on inverting M')
       return
    end if

    ! Calculate the partition sizes of the matrix problem
    n1 = no / 2
    n2 = no - n1

    if ( Npiv < max(n1,n2) ) then
       call die('Error in initialization of the pivoting &
            &array.')
    end if

    ! Point to the partitions
    sA1 = 1
    eA1 = n1 ** 2
    sB1 = eA1 + 1
    eB1 = sB1 + n1 * n2 - 1
    sC2 = eB1 + 1
    eC2 = sC2 + n1 * n2 - 1
    sA2 = eC2 + 1
    eA2 = sA2 + n2 ** 2 - 1
    A1 => zwork(sA1:eA1)
    B1 => zwork(sB1:eB1)
    C2 => zwork(sC2:eC2)
    A2 => zwork(sA2:eA2)

    ! Copy over the parts to ready the inversion
    idx = 0
    i1 = 0
    i2 = 0
    do i = 1 , n1
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             A1(i1) = M(idx)
          else
             i2 = i2 + 1
             B1(i2) = M(idx)
          end if
       end do
    end do
    i1 = 0
    i2 = 0
    do i = n1+1 , no
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             C2(i1) = M(idx)
          else
             i2 = i2 + 1
             A2(i2) = M(idx)
          end if
       end do
    end do
    
    ! Now we have created the matrices in the correct space
    ! Calculate both X1 and Y2
    
    ! Copy over A2 array
    t1 => M(sA2:eA2)
    t1(:) = A2(:)
    t2 => M(sB1:eB1)
    t2(:) = B1(:)
    ! Calculate X1/C2 (store in B1 in original matrix)
    call zgesv(n2,n1,t1,n2,ipiv,t2,n2,i)
    if ( i /= 0 ) call die('Error on inverting X1/C2')

    ! Copy over A1 array
    t1 => M(sA1:eA1)
    t1(:) = A1(:)
    t2 => M(sC2:eC2)
    t2(:) = C2(:)
    ! Calculate Y2/B1
    call zgesv(n1,n2,t1,n1,ipiv,t2,n1,i)
    if ( i /= 0 ) call die('Error on inverting Y2/B1')


! <<<<<<<<<<<<< direct method
    ! Calculate the diagonal inverted matrix
    ! Calculate: A1 - X1
    ! Copy over A1 array
!    t1 => M(sA1:eA1)
!    t1(:) = A1(:)
!    t2 => M(sB1:eB1) ! this is X1/C2
!    call zgemm('N','N',n1,n1,n2, &
!         zm1, C2,n1, t2,n2,z1, t1,n1)

!    call EYE(n1,A1)
!    call zgesv(n1,n1,t1,n1,ipiv,A1,n1,i)
!    if ( i /= 0 ) call die('Error on inverting A1-X1')
! >>>>>>>>>>>>>>>>>>>> end of direct method

    ! Calculate the diagonal inverted matrix
    ! Calculate: A1 - X1
    ! Copy over A1 array
    t2 => M(sB1:eB1) ! this is X1/C2
    call zgemm('N','N',n1,n1,n2, &
         zm1, C2,n1, t2,n2,z1, A1,n1)

    t1 => M(sA1:eA1)
    call mat_invert(A1,t1,n1)

! <<<<<<<<<<<<< direct method
    ! Calculate: A2 - Y2
    ! Copy over A2 array
!    t1 => M(sA2:eA2)
!    t1(:) = A2(:)
!    t2 => M(sC2:eC2) ! this is Y2/B1
!    call zgemm('N','N',n2,n2,n1, &
!         zm1, B1,n2, t2,n1,z1, t1,n2)

!    call EYE(n2,A2)
!    call zgesv(n2,n2,t1,n2,ipiv,A2,n2,i)
!    if ( i /= 0 ) call die('Error on inverting A2-Y2')
! >>>>>>>>>>>>>>>>>>>> end of direct method

    ! Calculate: A2 - Y2
    ! Copy over A2 array
    t2 => M(sC2:eC2) ! this is Y2/B1
    call zgemm('N','N',n2,n2,n1, &
         zm1, B1,n2, t2,n1,z1, A2,n2)

    t1 => M(sA2:eA2)
    call mat_invert(A2,t1,n2)

    ! Calculate the off-diagonal arrays
    ! Do matrix-multiplication
    ! Calculate: X1/C2 * M11
    t1 => M(sB1:eB1) ! this is X1/C2
    call zgemm('N','N',n2,n1,n1, &
         zm1, t1,n2, A1,n1,z0, B1,n2)

    ! Calculate the off-diagonal arrays
    ! Do matrix-multiplication
    ! Calculate: Y2/B1 * M22
    t1 => M(sC2:eC2) ! this is Y2/B1
    call zgemm('N','N',n1,n2,n2, &
         zm1, t1,n1, A2,n2,z0, C2,n1)

    ! Copy back the result
    ! Notice that the inverted matrix is back in the same matrix
    idx = 0
    i1 = 0
    i2 = 0
    do i = 1 , n1
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             M(idx) = A1(i1)
          else
             i2 = i2 + 1
             M(idx) = B1(i2)
          end if
       end do
    end do
    i1 = 0
    i2 = 0
    do i = n1+1 , no
       do j = 1 , no
          idx = idx + 1
          if ( j <= n1 ) then
             i1 = i1 + 1
             M(idx) = C2(i1)
          else
             i2 = i2 + 1
             M(idx) = A2(i2)
          end if
       end do
    end do

  end subroutine mat_invert


  ! We initialize the pivoting array for rotating the inversion
  subroutine init_mat_inversion(no)
    use alloc, only : re_alloc
    integer, intent(in) :: no
    Npiv = no / 2
    Npiv = max(no-Npiv,Npiv)
    Npiv = max(N_MAX,Npiv)

    ! Allocate space for the pivoting array
    call re_alloc(ipiv,1, Npiv, &
         name="mat_piv",routine='MatInversion')

  end subroutine init_mat_inversion

  subroutine clear_mat_inversion()
    use alloc, only: de_alloc
    Npiv = 0
    ! Deallocate the pivoting array
    call de_alloc(ipiv, &
         name="mat_piv",routine='MatInversion')
  end subroutine clear_mat_inversion

end module m_mat_invert
