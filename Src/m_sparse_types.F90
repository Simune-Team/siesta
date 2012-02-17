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
    character(len=256) :: name = "null_sparsity"
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
   character(len=256) :: name = "null_matrix"
   type(sparsity_t), pointer    :: sparsity => null()
   real(dp),         pointer    :: val(:,:) => null() ! Nonzero-element values
 end type

  ! For debugging purposes, to avoid memory leaks
  !integer, parameter                            :: nmax_sparsities = 5
  !type(sparsity_t), dimension(nmax_sparsities)  :: sparsity_pool

  !integer, save                                 :: n_sparsities = 0

CONTAINS

 subroutine new_sparsity(sparsity,nrows,ncols,nnzs,num,listptr,list,name)

   type(sparsity_t), pointer     :: sparsity
   integer, intent(in)  :: nrows, ncols, nnzs
   integer, intent(in)  :: num(:), listptr(:)
   integer, intent(in)  :: list(:)
   character(len=*), intent(in) :: name

   integer :: stat

   ! We release the previous incarnation
   ! This means that we relinquish access to the previous
   ! memory location. It will be deallocated when nobody
   ! else is using it.
   call release_sparsity(sparsity)
   !
   ! We allocate a new area for this incarnation
   allocate(sparsity,stat=stat)
   if (stat /= 0) then
      call die("Allocation error in new_sparsity")
   endif

   call re_alloc( sparsity%n_col, 1,nrows)
   call re_alloc( sparsity%list_ptr, 1,nrows)

   sparsity%nrows = nrows
   sparsity%ncols = ncols
   sparsity%nnzs  = nnzs
   sparsity%n_col(1:nrows) = num(1:nrows)
   sparsity%list_ptr(1:nrows) = listptr(1:nrows)

   if (nnzs /= sum(num(1:nrows))) then
      call die("nnzs mismatch in new_sparsity")
   endif

   call re_alloc( sparsity%list_col, 1,nnzs)
   sparsity%list_col(1:nnzs) = list(1:nnzs)

   sparsity%refcount = 0
   sparsity%initialized = .true.   
   sparsity%name = trim(name)

   print *, "--> allocated sparsity: " // trim(sparsity%name)
   ! Should we stake a claim to this record?
   ! Not if the sparsity is only going to be used in matrices...
   ! call retain_sparsity(sparsity)
   
 end subroutine new_sparsity

 subroutine release_sparsity(sparsity)
   type(sparsity_t), pointer     :: sparsity

   if (associated(sparsity)) then
      sparsity%refcount = sparsity%refcount - 1
      print *, "--> released sparsity: " //  &
                trim(sparsity%name) //", new refcount: ", sparsity%refcount
      if (sparsity%refcount == 0) then
         call de_alloc( sparsity%n_col)
         call de_alloc( sparsity%list_ptr)
         call de_alloc( sparsity%list_col)
         print *, "--> deallocated sparsity: " // trim(sparsity%name) 
         deallocate(sparsity)
         nullify(sparsity)
      endif
   endif
   
 end subroutine release_sparsity

 subroutine retain_sparsity(sparsity)
   type(sparsity_t), pointer     :: sparsity

   if (associated(sparsity)) then
      sparsity%refcount = sparsity%refcount + 1
      print *, "--> retained sparsity: " //  &
                trim(sparsity%name) //", new refcount: ", sparsity%refcount
   else
      call die("Attempted to retain unset sparsity")
   endif
   
 end subroutine retain_sparsity

 function matching_sparsity(os,ns) result (ok)
   type(sparsity_t), pointer :: os, ns
   logical                   :: ok

    ok = associated(os,ns)
  end function matching_sparsity

!-----------------------------------------------------------------------
 subroutine new_sparse_matrix( matrix, sparsity, dim2, name, fill_value )

   implicit none

   type(sparseMatrix_t), intent(out)      :: matrix
   type(sparsity_t), pointer              :: sparsity 
   integer, intent(in)                    :: dim2  ! Should it be lbound/ubound?
   character(len=*), intent(in)           :: name
   real(dp), intent(in), OPTIONAL         :: fill_value


   if (.not. associated(sparsity)) then
      call die( "sparsity not initialized in new_sparse_matrix")
   endif
   !
   call release_sparse_matrix(matrix)   ! In case it was holding info
   !
   matrix%sparsity => sparsity
   matrix%name = trim(name)
   print *, "Allocating info for matrix " // trim(matrix%name) // &
               " with sparsity " // trim(matrix%sparsity%name)
   call retain_sparsity(matrix%sparsity)
   call re_alloc( matrix%val, 1,matrix%sparsity%nnzs, 1, dim2, &
                  routine='new_sparse_matrix', &
                  shrink=.true., copy=.false.)
   !
   if (present(fill_value)) then
      matrix%val(:,:) = fill_value   ! Optional
   endif

 end subroutine new_sparse_matrix

!--------------------------------------------------------------------------
 subroutine copy_sparse_matrix (destination, name, source)
   type(sparseMatrix_t),     intent(in)     :: source
   type(sparseMatrix_t),     intent(inout)  :: destination
   character(len=*), intent(in)             :: name

   type(sparsity_t), pointer :: os, ns
   
   if (associated(source%sparsity)) then
      os => source%sparsity
   else
      call die("Source matrix " // trim(source%name) // &
               " not associated in copy_sparse_matrix")
   endif

   print *, "Starting copy of matrix " // trim(source%name) // &
               " into new matrix " // trim(name)

   if (associated(destination%sparsity)) then
      ! Destination matrix already setup
      ns => destination%sparsity
      if (matching_sparsity(os,ns)) then
         ! Both matrices are compatible
         ! simply copy the data, but check first
         if (size(source%val,dim=2) /= size(destination%val,dim=2)) then
            ! This case has to be dealt with 
            ! by explicit mapping of the data
            call die("dim2 mismatch in copy")
         else
            destination%val(:,:) = source%val(:,:)
         endif
         print *, "Copied data into existing matrix " // trim(destination%name) // &
               ". New name: " // trim(name)
         destination%name = trim(name)
      else
         call new_sparse_matrix(destination,source%sparsity, &
                                size(source%val,dim=2),name)
         destination%val(:,:) = source%val(:,:)

      endif

   else  ! We need to set up the new matrix

      call new_sparse_matrix(destination,source%sparsity, &
                             size(source%val,dim=2),name)
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

   ! Print messages only if needed
   if (associated(matrix%sparsity)) then
      print *, "Releasing  matrix " // trim(matrix%name) 
      if (associated(matrix%val)) then
         call de_alloc( matrix%val )
      else
         call die("Data for matrix was not allocated...")
      endif
      matrix%name = "null_matrix"
   endif
   call release_sparsity(matrix%sparsity)

 end subroutine release_sparse_matrix

END MODULE m_sparse_types
