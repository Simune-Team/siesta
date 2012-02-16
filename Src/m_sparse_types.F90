MODULE m_sparse_types

 use precision, only: dp
 use sys,       only: die
 use alloc,     only: re_alloc, de_alloc

 implicit none

PUBLIC :: sparsity_t, sparseMatrix_t
PUBLIC :: new_sparsity
PUBLIC :: new_sparse_matrix, copy_sparse_matrix, release_sparse_matrix
PUBLIC :: map_indexes_sparse_matrix, map_data_sparse_matrix

PRIVATE ! Nothing is declared public beyond this point

 type sparsity_t
    private
!   integer        :: n_row_g=0  ! Global number of rows
!   integer        :: n_col_g=0  ! Global number of columns
    integer        :: nrows
    integer        :: ncols
    integer        :: nnzs
!   integer,pointer:: n_row_node(:)=>null() ! # rows in each node
!   integer,pointer:: row_g2l(:)   =>null() ! Global to local row index
!   integer,pointer:: row_l2g(:)   =>null() ! Local to global row index
   integer,pointer:: n_col(:)     =>null() ! Nonzero cols of each row
   integer,pointer:: list_col(:)  =>null() ! Index of nonzero columns
   integer,pointer:: list_ptr(:)  =>null() ! First element of each row
   integer        :: refcount = 0          ! Reference count
   logical        :: initialized = .false.
 end type

 type sparseMatrix_t
    private
   type(sparsity_t), pointer    :: sparsity => null()
   real(dp),         pointer    :: val(:,:) => null() ! Nonzero-element values
 end type

CONTAINS

 subroutine new_sparsity(sparsity,nrows,ncols,nnzs,num,listptr,list)

   type(sparsity_t), pointer     :: sparsity
   integer, intent(in)  :: nrows, ncols, nnzs
   integer, intent(in)  :: num(:), listptr(:)
   integer, intent(in)  :: list(:)

   integer :: stat

   if (associated(sparsity)) then
      call die( "sparsity already allocated")
   endif
   allocate(sparsity,stat=stat)
   call re_alloc( sparsity%n_col, 1,nrows)
   call re_alloc( sparsity%list_ptr, 1,nrows)

   sparsity%nrows = nrows
   sparsity%ncols = ncols
   sparsity%nnzs  = nnzs
   sparsity%n_col(1:nrows) = num(1:nrows)
   sparsity%list_ptr(1:nrows) = listptr(1:nrows)

   if (nnzs /= sum(num(1:nrows))) then
      call die("nnzs mismatch")
   endif

   call re_alloc( sparsity%list_col, 1,nnzs)
   sparsity%list_col(1:nnzs) = list(1:nnzs)

   sparsity%refcount = 0
   sparsity%initialized = .true.
   
 end subroutine new_sparsity

 subroutine release_sparsity(sparsity)
   type(sparsity_t), pointer     :: sparsity

   if (associated(sparsity)) then
      sparsity%refcount = sparsity%refcount - 1
      if (sparsity%refcount == 0) then
         call de_alloc( sparsity%n_col)
         call de_alloc( sparsity%list_ptr)
         call de_alloc( sparsity%list_col)
         deallocate(sparsity)
         nullify(sparsity)
      endif
   endif
   
 end subroutine release_sparsity

 subroutine retain_sparsity(sparsity)
   type(sparsity_t), pointer     :: sparsity

   if (associated(sparsity)) then
      sparsity%refcount = sparsity%refcount + 1
   else
      STOP "sparsity not retained"
   endif
   
 end subroutine retain_sparsity

 function matching_sparsity(os,ns) result (ok)
   type(sparsity_t), pointer :: os, ns
   logical                   :: ok

    ok = associated(os,ns)
  end function matching_sparsity

!-----------------------------------------------------------------------
 subroutine new_sparse_matrix( matrix, sparsity, dim2, value )

   implicit none

   type(sparseMatrix_t), intent(out)      :: matrix
   type(sparsity_t), pointer              :: sparsity 
   integer, intent(in)                    :: dim2  ! Should it be lbound/ubound?
   real(dp), intent(in), OPTIONAL         :: value


   if (.not. associated(sparsity))  call die( "sparsity not initialized")
   call re_alloc( matrix%val, 1,sparsity%nnzs, 1, dim2, routine='sparse_alloc' )
   if (present(value)) then
      matrix%val(:,:) = value   ! Optional
   endif
   matrix%sparsity => sparsity
   call retain_sparsity(matrix%sparsity)

 end subroutine new_sparse_matrix

!--------------------------------------------------------------------------
 subroutine copy_sparse_matrix (destination, source)
   type(sparseMatrix_t),     intent(in)     :: source
   type(sparseMatrix_t),     intent(inout)  :: destination

   type(sparsity_t), pointer :: os, ns
   
   if (associated(source%sparsity)) then
      os => source%sparsity
   else
      STOP "source matrix not associated"
   endif

   if (associated(destination%sparsity)) then
      ! Destination matrix already setup
      ns => destination%sparsity
      if (matching_sparsity(os,ns)) then
         ! Both matrices are compatible
         ! simply copy the data, but check first
         if (size(source%val,dim=2) /= size(destination%val,dim=2)) then
            ! We will deal with this case later
            STOP "dim2 mismatch in copy"
         else
            destination%val(:,:) = source%val(:,:)
         endif
      else
         call release_sparsity(ns)
         destination%sparsity=>os
         call retain_sparsity(destination%sparsity)
         call de_alloc(destination%val)
         call re_alloc(destination%val,1,os%nnzs,1,size(source%val,dim=2),routine="..")
         destination%val(:,:) = source%val(:,:)
      endif

   else  ! We need to set up the new matrix

      call new_sparse_matrix(destination,source%sparsity,size(source%val,dim=2))
      destination%val(:,:) = source%val(:,:)

   endif

 end subroutine copy_sparse_matrix

 subroutine map_data_sparse_matrix(matrix,ptr)
   type(sparseMatrix_t), intent(in)    :: matrix
   real(dp), pointer                   :: ptr(:,:)

   if (.not.associated(matrix%sparsity)) then
      call die("Matrix not associated -- cannot map data")
   else
      ptr => matrix%val
   endif
 end subroutine map_data_sparse_matrix

 subroutine map_indexes_sparse_matrix(matrix,nrows,ncols,nnzs,num,listptr,list)
   type(sparseMatrix_t), intent(in)    :: matrix
   integer, intent(out), OPTIONAL      :: nrows, ncols, nnzs
   integer, pointer, OPTIONAL          :: num(:), listptr(:), list(:)

   type(sparsity_t), pointer :: sp

   if (.not.associated(matrix%sparsity)) then
      call die("Matrix not associated -- cannot map data")
   else
      sp => matrix%sparsity
      if (present(nrows)) then
         nrows = sp%nrows
      endif
      if (present(ncols)) then
         ncols = sp%ncols
      endif
      if (present(nnzs)) then
         nnzs = sp%nnzs
      endif
      if (present(num)) then
         num   => sp%n_col
      endif
      if (present(listptr)) then
         listptr => sp%list_ptr
      endif
      if (present(list)) then
         list    => sp%list_col
      endif
   endif

 end subroutine map_indexes_sparse_matrix

!-----------------------------------------------------------------------
 subroutine release_sparse_matrix( matrix )

   implicit none

   type(sparseMatrix_t), intent(inout)      :: matrix

   if (.not. associated(matrix%sparsity)) then
      call die( "matrix not initialized -- cannot release")
   endif
   call release_sparsity(matrix%sparsity)
   call de_alloc( matrix%val )

 end subroutine release_sparse_matrix

END MODULE m_sparse_types
