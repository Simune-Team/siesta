subroutine Inver( a, b, n )

! Subroutine for matrix inversion
! Contributed by William Mattson (wmattson@uiuc.edu), April 2001
! 
! Input :
!
!  integer          n      : Dimension of matrices
!  double precision a(n,n) : Matrix to be inverted
!
! Output :
!
!  double precision b(n,n) : Inverse of a
!

  implicit none

  integer, intent( in ) :: n
  double precision, intent( in ) :: a( n, n )
  double precision, intent( out ) :: b( n, n )
  double precision :: biggest, invPivot, scale
  double precision, dimension( n ) :: temp
  integer, dimension( n ) :: pivotAvailable, rowIndex, colIndex
  integer :: i, j, k, row, col, nPivotAvailable
  double precision, parameter :: zero = 0.0d0, one = 1.0d0

  do i = 1, n
    pivotAvailable( i ) = i
  end do
  b = a
  i = 0

  ! Loop over all rows for pivoting.
  do nPivotAvailable = n, 1, -1

    ! Find the row and column for the element with the biggest value among
    ! the rows and columns that haven't already been used as pivots.
    i = i + 1
    biggest = zero
    do j = 1, nPivotAvailable
      do k = 1, nPivotAvailable
        if( abs( b( pivotAvailable(j), pivotAvailable(k) )) .ge. biggest ) then
          biggest = abs( b( pivotAvailable(j), pivotAvailable(k) ))
          row = pivotAvailable(j)
          col = pivotAvailable(k)
        end if
      end do
    end do

    ! If the biggest value left is zero the matrix is singular.
    if( biggest .eq. zero ) stop 'Singular Matrix in Inver'

    ! Remove the current pivot column from the
    ! list of available rows and columns
    pivotAvailable( col:(n-1) ) = pivotAvailable( (col+1):n )

    ! Swap rows to get the biggest value on the diagonal.
    temp = b( row, : )
    b( row, : ) = b( col, : )
    b( col, : ) = temp

    ! Store the information about the original location of the
    ! row and column so that we can unswap them later
    rowIndex( i ) = row
    colIndex( i ) = col

    ! Invert current row.
    invPivot = one / b( col, col )
    b( col, col ) = one
    b( col, : ) = b( col, : ) * invPivot

    ! Adjust the other rows.
    do row = 1, n
      if( row .ne. col ) then
        scale = b( row, col )
        b( row, col ) = zero
        b( row, : ) = b( row, : ) - b( col, : ) * scale
      end if
    end do
  end do

  ! Unswap the columns in the matrix that correspond
  ! to the rows we swapped earlier.
  do col = n, 1, -1
    if( rowIndex( col ) .ne. colIndex( col )) then
      temp = b( :, rowIndex( col ))
      b( :, rowIndex( col )) = b( :, colIndex( col ))
      b( :, colIndex( col )) = temp
    end if
  end do

end subroutine Inver
