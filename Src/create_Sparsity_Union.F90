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

! This module was intended for the use of the transiesta
! package. It allows to alter a sparse array by adding a 
! dense region to an existing sparsity pattern

! Generally for Transiesta you need the full electrodes
! to contain the self-energies.
! this can be introduced by something like this:

!    call crtSparsity_Union(block_dist,ts_sp,&
!         no_BufL+1, no_BufL+1, &
!         no_L, no_L, & ! insertion of the block
!         sp_uc)
!    if(IONode)call print_type(sp_uc)

!    call crtSparsity_Union(block_dist,sp_uc,&
!         no_BufL+no_u_LCR-no_R+1, no_BufL+no_u_LCR-no_R+1, &
!         no_R, no_R, & ! insertion of the block
!         tsup_sp)
!    if(IONode)call print_type(tsup_sp)

module create_Sparsity_Union
  
  use precision, only : dp

  implicit none

  private

  public :: crtSparsity_Union

contains

  ! Subroutine will create a new pattern which matches will be the union of 
  ! the 'in' sparse pattern and a segmented block in the full matrix.
  !
  ! Then routine name reflects its use
  ! crt for "create"
  ! Sparsity for "Sparsity" and
  ! Union for uniting the sparsity patterns.
  ! 
  ! Needless to say, this routine makes not much sense when the final sparsity is > 50%
  ! Also you should consider that this will only work properly on unit-cell sparse
  ! patterns.
  subroutine crtSparsity_Union(dit,in,in_i,in_j,ni,nj,out)
    use class_Sparsity
    use class_OrbitalDistribution
    use parallel, only : Node

    ! The sparse distribution
    type(OrbitalDistribution), intent(inout) :: dit
    ! The matrix which contains the unit cell.
    type(Sparsity), intent(in out) :: in
    ! The row, col, # of rows, # of columns for the dense part.
    integer, intent(in) :: in_i, in_j, ni, nj
    ! Maybe (in) is not needed, however, for this, it does't make a
    ! difference, everything is overwritten
    type(Sparsity), intent(in out) :: out
    
    ! We need space to create the new sparsity pattern:
    integer :: n_rows, n_rows_g, n_nzs
    integer, allocatable :: num(:), listptr(:), list(:)

    integer :: lir,ir,i,iu,ptr,l,j

    ! Local variables...
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)
    logical :: dense_inserted

    ! Save the rows ( this is the same for all cases)
    ! Even for TM which typically have 0 entries
    ! in some rows. However, it provides the full information.
    n_rows   =  nrows   (in)
    n_rows_g =  nrows_g (in)
    ncol     => n_col   (in)
    l_ptr    => list_ptr(in)
    l_col    => list_col(in)

    ! We no a check of arguments before anything...
    if ( in_i < 1 .or. n_rows_g < in_i + ni - 1 ) then
       write(*,'(a,2(i0,tr1),a,i0)') 'The rows requested is not within &
            &the sparse pattern: ',in_i,in_i+ni-1,'vs. ',n_rows_g
       call die('Unifying a sparse matrix and a dense requires the &
            &dense part to be a subset of the sparse matrix. &
            &This is not enforced.')
    end if
    if ( in_j < 1 .or. n_rows_g < in_j + nj - 1 ) then
       write(*,'(a,2(i0,tr1),a,i0)') 'The columns requested is not within &
            &the sparse pattern: ',in_j,in_j+nj-1,'vs. ',n_rows_g
       call die('Unifying a sparse matrix and a dense requires the &
            &dense part to be a subset of the sparse matrix. &
            &This is not enforced.')
    end if

    ! Prepare creation of num and listptr arrays
    allocate(num    (n_rows))
    allocate(listptr(n_rows))

    ! The list pointer for the first entry is always the "first"
    ! element, hence we can already initialize it here...
    listptr(1) = 0

    ! We initialize the sparsity creation via the first row
    do lir = 1 , n_rows

       ! Retrieve the global row-index
       ir = index_local_to_global(dit,lir,Node)

       ! The easy part is if we are not in the row containing the dense part
       if ( ir < in_i .or. in_i + ni <= ir ) then
          num(lir) = ncol(lir)
       else
          ! We assume a block-cyclic distribution
          ! Which means that the current row has the dense part
          ! We know all the dense elements, then we only need to
          ! count the non-dense
          num(lir) = nj
          ! Retrieve the pointer to the original sparsity
          ptr = l_ptr(lir)
          do i = 1 , ncol(lir)
             ! we check whether the sparse element is
             ! outside the dense part (in that case we need
             ! to count it)
             if ( l_col(ptr+i) <  in_j           .or. &
                  in_j + nj    <= l_col(ptr+i) ) then
                num(lir) = num(lir) + 1
             end if
          end do
       end if

       ! Update list pointer
       if ( lir > 1 ) &
            listptr(lir) = listptr(lir-1) + num(lir-1)
    end do

    ! The number of non-zero elements in the requested array
    n_nzs = listptr(n_rows) + num(n_rows)

    ! Now we will create the list array...
    allocate(list(n_nzs))

    do lir = 1, n_rows

       ! Retrieve the current pointer
       iu = listptr(lir)

       ! The easy part is if num(lir) == n_col(in,lir)
       ! in which case, we can be sure that we are outside the dense
       ! region
       if ( num(lir) == ncol(lir) ) then
          ! Copy over values from the original sparse pattern
          ptr = l_ptr(lir)
          do i = 1 , num(lir)
             list(iu+i) = l_col(ptr+i)
          end do
       else
          ! Initialize original sparse pattern pointer
          ptr = l_ptr(lir)
          ! We count the number of entries and compare in the end
          ! this way we ensure that we at least count the 
          ! correct number of entries
          j = 0
          ! We only need to insert the dense matrix once
          ! (this should release the condition of a non-sorted
          !  list_col)
          dense_inserted = .false.
          ! We loop over the local sparse pattern
          do i = 1 , ncol(lir)
             ! When we are on either side of the dense part
             ! we simply copy the column index
             if ( l_col(ptr+i) < in_j .or. &
                  nj+in_j <= l_col(ptr+i) ) then
                j = j + 1
                list(iu+j) = l_col(ptr+i)
             else if ( .not. dense_inserted ) then

                ! we are within the dense part of the array
                dense_inserted = .true.
                ! insert the dense format
                do l = 0 , nj - 1
                   j = j + 1
                   list(iu+j) = in_j + l
                end do

             end if

          end do
          
          ! In the rare occasion that no elements
          ! exists in the dense region, we pad the
          ! back of the list arry with the dense matrix elements
          if ( .not. dense_inserted ) then
             ! insert the dense format
             do l = 0 , nj - 1
                j = j + 1
                list(iu+j) = in_j + l
             end do
          end if
          
          if ( j /= num(lir) ) then
             call die('Creating the column indices was &
                  &ill-formatted. Please ask the developers for help.')
          end if
       end if
    end do

    call newSparsity(out,n_rows,n_rows_g,n_nzs,num,listptr,list, &
         name='(Dense-Union region of: '//name(in)//')', &
         ncols=ncols(in),ncols_g=ncols_g(in))

    ! We need to deallocate the arrays used, remember, that they
    ! are allocated.
    ! The newSparsity copies the values...
    deallocate(num,listptr,list)

  end subroutine crtSparsity_Union

end module create_Sparsity_Union
