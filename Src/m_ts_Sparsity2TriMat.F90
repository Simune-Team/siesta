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
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.
! * It has been heavily inspired by the original authors of the 
!   Transiesta code (hence the references here are still remaining) *

! Module for converting a sparsity pattern to an
! "optimal" tri-diagonal matrix...

! Routine for converting a SIESTA sparsity pattern to a tri-diagonal form
! It requires that the sparsity pattern is fully contained in the current
! processor.

! This module has been fully developed by Nick Papior Andersen, 2013
! nickpapior@gmail.com

! Please contact the author before utilization in other routines.

module m_ts_Sparsity2TriMat
  
  implicit none

  private

  public :: ts_Sparsity2TriMat

  integer, parameter :: VALID = 0
  integer, parameter :: NONVALID_SIZE = 1
  integer, parameter :: NONVALID_ELEMENT_CONTAIN = 2
  integer, parameter :: NONVALID_TS_ELECTRODE = 3
  
contains

  ! IF parts == 0 will create new partition
  subroutine ts_Sparsity2TriMat(sp,parts,n_part)
    use class_Sparsity
    use parallel, only : IONode
    use alloc, only : re_alloc, de_alloc
    use m_ts_electype
    use m_ts_options, only : no_BufL, no_BufR
    use m_ts_options, only : ElLeft, ElRight

    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! The sizes of the parts in the tri-diagonal matrix
    integer, intent(in out) :: parts
    integer, pointer :: n_part(:)
    ! Local variables
    integer, pointer :: tmp_part(:)
    integer :: i, N, val

    ! Establish a guess on the partition of the tri-diagonal 
    ! matrix...
    if ( parts == 0 ) then
       parts = 0
       N = 0
       do while ( N < nrows_g(sp) - no_BufL - no_BufR)
          ! Albeit this is "slow" we should never exceed 1000 re-allocations
          ! (that would mean a huge system)
          parts = parts + 1
          call re_alloc(n_part, 1, parts, copy=.true.)
          if ( parts == 1 ) then
             call guess_end_part(sp, parts, n_part,first=.true.)
          else
             call guess_next_part_size(sp, parts, parts, n_part)
          end if
          N = N + n_part(parts)
       end do
    end if

    if ( parts < 3 ) then
       if ( IONode ) then 
          write(*,'(a)') 'Could not determine an optimal tri-diagonalization &
               &partition'
          write(*,'(a,i0)') 'Found: ',parts
          write(*,'(1000000(tr1,i0))') n_part
       end if
       parts = 3 
       call re_alloc(n_part, 1, parts)
       n_part(1) = TotUsedOrbs(ElLeft)
       n_part(3) = TotUsedOrbs(ElRight)
       n_part(2) = nrows_g(sp) - no_BufL - no_BufR - n_part(1) - n_part(3)

    else

       ! Even out the matrix part sizes...

       call re_alloc(tmp_part,1,parts)
       N = 1
       do while ( N /= 0 )
          tmp_part(:) = n_part(:)
          do i = 1 , parts
             call even_out_parts(sp, parts, n_part, i)
          end do
          N = maxval(abs(tmp_part-n_part))
       end do
       call de_alloc(tmp_part)

    end if

    ! The parts now have almost the same size and we will check that it
    ! is a valid thing, if not, we will revert to the other method of
    ! creating the tri-diagonal sparsity pattern
    ! We do not expect this to fail. After all we check that we can
    ! even out the partitions.
    ! The most probable thing is that the electrodes are not
    ! contained in the first two parts.
    if ( ts_valid_tri(sp, parts, n_part) /= VALID ) then
       write(*,'(a,i0)') 'Current parts: ',parts
       write(*,'(10000000(tr1,i0))') n_part
       call die('Contact the developers. (missing implementation). &
            &You appear to have a special form of electrode.')
    end if

  end subroutine ts_Sparsity2TriMat

  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in n_part(part)
  subroutine guess_next_part_size(sp,part,parts,n_part)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the part we are going to create
    integer, intent(in) :: part, parts
    integer, intent(inout) :: n_part(parts)
    ! Local variables
    integer :: i, sRow, eRow, mcol
    
    ! We are now checking a future part
    ! Hence we must ensure that the size is what
    ! is up to the last parts size, and thats it...
    sRow = 1
    do i = 1 , part - 2
       sRow = sRow + n_part(i)
    end do
    eRow = sRow + n_part(part-1) - 1
    ! We will check in between the above selected rows and find the 
    ! difference in size...
    n_part(part) = 0
    do i = sRow, eRow
       ! this is the # of elements from the RHS of the 'part-1'
       ! part of the tridiagonal matrix and out to the last element of
       ! this row...
       mcol = max_col(sp,i) - eRow
       if ( n_part(part) < mcol ) then
          n_part(part) = mcol
       end if
    end do

  end subroutine guess_next_part_size

  subroutine guess_previous_part_size(sp,part,parts,n_part)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the part we are going to create
    integer, intent(in) :: part, parts
    integer, intent(inout) :: n_part(parts)
    ! Local variables
    integer :: i, sRow, eRow, mcol, ncol
    
    ! We are now checking a future part
    ! Hence we must ensure that the size is what
    ! is up to the next parts size, and thats it...
    ncol = ncols(sp) + 1
    eRow = nrows_g(sp)
    do i = parts , part + 2 , -1
       eRow = eRow - n_part(i)
    end do
    sRow = eRow - n_part(part+1) + 1

    ! We will check in between the above selected rows and find the 
    ! difference in size...
    n_part(part) = 0
    do i = sRow, eRow
       ! this is the # of elements from the LHS of the 'part+1'
       ! part of the tridiagonal matrix and out to the firs element of
       ! this row...
       mcol = sRow - min_col(sp,i)
       if ( n_part(part) < mcol ) then
          n_part(part) = mcol
       end if
    end do
    
  end subroutine guess_previous_part_size

  ! We will guess the size of this (part) tri-diagonal part by
  ! searching for the size of the matrix that matches that of the previous 
  ! part.
  ! We save it in n_part(part)
  subroutine guess_end_part(sp,parts,n_part,first)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the part we are going to create
    integer, intent(in) :: parts
    integer, intent(inout) :: n_part(parts)
    logical, intent(in) :: first
    ! Local variables
    integer :: i, mcol, ncol
    
    ! When we need to guess the first part size, then
    ! we will check the first row and select half the 
    ! maximum column. This makes sense as then the tri-diagonal
    ! parts would be the same size throughout the matrix.

    ! Then we will check each of those columns and
    ! adjust the size accordingly...
    if ( first ) then
       mcol = max_col(sp,1)
       
       ! We now have an initial guess of the number 
       ! of rows in the first part
       ! Check that it "makes sense".
       n_part(1) = mcol / 2 + mod(mcol,2)
       i = 2
       do while ( i <= n_part(1) ) 
          if ( max_col(sp,i) > mcol ) then
             mcol = max_col(sp,i)
             n_part(1) = mcol / 2 + mod(mcol,2)
          end if
          i = i + 1
       end do

    else
       ncol = ncols(sp) + 1
       mcol = ncol - min_col(sp,nrows(sp))

       ! We now have an initial guess of the number 
       ! of rows in the first part
       ! Check that it "makes sense".
       n_part(parts) = mcol / 2 + 1

       i = nrows(sp) - 1
       do while ( nrows(sp) - i <= n_part(parts) )
          if ( ncol - min_col(sp,i) > mcol ) then
             mcol = ncol - min_col(sp,i)
             n_part(parts) = mcol / 2 + 1
          end if
          i = i - 1
       end do
       
    end if

    ! We now have guessed the first/last part
  end subroutine guess_end_part

  subroutine even_out_parts(sp,parts,n_part, n)
    use class_Sparsity
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the part we are going to create
    integer, intent(in) :: parts
    integer, intent(in out) :: n_part(parts)
    integer, intent(in) :: n
    ! Local variables
    integer :: copy_n_part
    integer :: sRow, eRow
    integer :: i

    if ( parts < 2 ) call die('You cannot use tri-diagonalization &
         &without having at least 2 parts')
    
    copy_n_part = n_part(n) + 1

    sRow = 1
    do i = 1 , n - 1
       sRow = sRow + n_part(i)
    end do
    eRow = sRow + n_part(n) - 1

    ! We will continue to shift columns around
    ! until we do not shift columns any more...
    do while ( n_part(n) - copy_n_part /= 0 )

       ! Copy the current partition so that we can check in the
       ! next iteration...
       copy_n_part = n_part(n)

       ! TODO - consider adding a flag for memory reduced utilization of TRI
       ! this will require the two electrodes parts to be larger
       ! than the other parts...

       if ( n == 1 ) then
          ! If we have the first part we can always shrink it
          call even_if_larger(sRow,n_part(n),n_part(n+1))
       else if ( n == parts ) then
          ! If we have the last part we can always shrink it
          call even_if_larger(eRow,n_part(n),n_part(n-1))
       else
          ! all middle parts have some requirements:
          
          ! 1. if you wish to shrink it left, then:
          !    the first row must not have any elements
          !    extending into the right part
          if ( max_col(sp,sRow) <= eRow ) then
             call even_if_larger(sRow,n_part(n),n_part(n-1))
          end if

          ! 2. if you wish to shrink it right, then:
          !    the last row must not have any elements
          !    extending into the left part
          if ( sRow <= min_col(sp,eRow) ) then
             call even_if_larger(eRow,n_part(n),n_part(n+1))
          end if

       end if

    end do

  contains

    subroutine even_if_larger(Row,p1,p2)
      integer , intent(inout) :: Row, p1, p2
      if ( p1 > p2 ) then
         p1 = p1 - 1
         p2 = p2 + 1
         Row = Row + 1
      end if
    end subroutine even_if_larger

  end subroutine even_out_parts

  function max_col(sp,row)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the row which we will check for
    integer, intent(in) :: row
    ! The result
    integer :: max_col, ptr, nr
    integer, pointer :: l_col(:)
    call retrieve(sp,list_col=l_col,nrows_g=nr)
    ptr     =  list_ptr(sp,row)
    max_col =  maxval(UCORB( &
         l_col(ptr+1:ptr+n_col(sp,row)),nr))
  end function max_col

  function min_col(sp,row)
    use class_Sparsity
    use geom_helper, only : UCORB
    ! The sparsity pattern
    type(Sparsity), intent(inout) :: sp
    ! the row which we will check for
    integer, intent(in) :: row
    ! The result
    integer :: min_col, ptr, nr
    integer, pointer :: l_col(:)
    call retrieve(sp,list_col=l_col,nrows_g=nr)
    ptr     =  list_ptr(sp,row)
    min_col =  minval(UCORB( &
         l_col(ptr+1:ptr+n_col(sp,row)),nr))
  end function min_col

  function valid_tri(sp,parts,n_part) result(val) 
    use class_Sparsity
    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: parts, n_part(parts)
    integer :: val
    ! Local variables
    integer :: i, ir, N, Nm1, Np1

    val = VALID
    ! Calculate the size of the tri-matrix
    N = sum(n_part)

    ! Easy check, if the number of rows
    ! does not sum up to the total number of rows.
    ! Then it must be invalid...
    if ( N /= nrows(sp) ) then
       val = NONVALID_SIZE
       return
    end if

    ! Check that every element is contained in the 
    ! tri-diagonal matrix...
    N = 1
    Nm1 = 1
    Np1 = n_part(1)
    do i = 1 , parts - 1
       ! Update the size of the part after this
       Np1 = Np1 + n_part(i+1)
       
       do ir = N , N + n_part(i) - 1
          if ( Nm1 > min_col(sp,ir) .or. &
               max_col(sp,ir) > Np1 ) then
             val = - ir 
             return
          end if
       end do
       ! Update loop
       N = N + n_part(i)

       if ( i > 1 ) then
          ! Update the previous part
          Nm1 = Nm1 + n_part(i-1)
       end if
    end do
       
  end function valid_tri

  ! Validation routine for the tri-diagonal splitting
  ! with a Transiesta tri-diagonal matrix.
  ! It will first check for the electrode size
  function ts_valid_tri(sp,parts,n_part) result(val)
    use class_Sparsity
    use m_ts_electype
    use m_ts_options, only: ElLeft, ElRight
    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: parts, n_part(parts)
    integer :: val

    if ( sum(n_part(1:2)) < TotUsedOrbs(ElLeft) .or. &
         sum(n_part(parts-1:parts)) < TotUsedOrbs(ElRight) ) then
       val = NONVALID_TS_ELECTRODE
       return
    end if

    val = valid_tri(sp,parts,n_part)
    if ( val /= VALID ) return
    
  end function ts_valid_tri


  subroutine ts_adjust_tri_size(sp,parts,n_part)
    use class_Sparsity
    use m_ts_electype
    use m_ts_options, only: ElLeft, ElRight
    type(Sparsity), intent(inout) :: sp
    integer, intent(in) :: parts
    integer, intent(inout) :: n_part(parts)
    integer :: no_L, no_R, no_C
    no_L = TotUsedOrbs(ElLeft)
    no_R = TotUsedOrbs(ElRight)
    no_C = sum(n_part) - no_L - no_R

    ! If there is only three parts then
    ! the routines are optimized for 
    ! electrode-central
    if ( parts == 3 ) then
       n_part(1) = no_L
       n_part(2) = no_C
       n_part(3) = no_R
       return
    end if

    do while ( no_L > sum(n_part(1:2)) )
       if ( n_part(1) < n_part(2) ) then
          n_part(1) = n_part(1) + 1
          n_part(3) = n_part(3) - 1
       else if ( n_part(1) > n_part(2) ) then
          n_part(2) = n_part(2) + 1
          n_part(3) = n_part(3) - 1
       else
          n_part(1) = n_part(1) + 1
          n_part(2) = n_part(2) + 1
          n_part(3) = n_part(3) - 2
       end if
    end do
    do while ( no_R > sum(n_part(parts-1:parts)) )
       if ( n_part(parts-1) < n_part(parts) ) then
          n_part(parts-1) = n_part(parts-1) + 1
          n_part(parts-2) = n_part(parts-2) - 1
       else if ( n_part(parts-1) > n_part(parts) ) then
          n_part(parts) = n_part(parts) + 1
          n_part(parts-2) = n_part(parts-2) - 1
       else
          n_part(parts-1) = n_part(parts-1) + 1
          n_part(parts) = n_part(parts) + 1
          n_part(parts-2) = n_part(parts-2) - 2
       end if
    end do
  end subroutine ts_adjust_tri_size

end module m_ts_Sparsity2TriMat

