MODULE m_sparse_types

! Derived types for sparse matrices and their sparsity structure
! J.M.Soler. Jan.2007
!
! Usage:
!   type(sparsity_t),target :: sparsity
!   type(sparseMatrix)      :: H, S
!   call set_my_sparsity( sparsity )
!   call sparse_alloc( H, sparsity )
!   call sparse_alloc( S, sparsity )

 use precision, only: dp
 use sys,       only: die
 use alloc,     only: re_alloc, de_alloc

 implicit none

PUBLIC :: sparsity_t, sparseMatrix, sparse_alloc
PUBLIC :: build_sparsity, change_sparsity, setup_old_style_arrays
PUBLIC :: copy_sparse_matrix, copy_sparsity

PRIVATE ! Nothing is declared public beyond this point

 type sparsity_t
!   integer        :: n_row_g=0  ! Global number of rows
!   integer        :: n_col_g=0  ! Global number of columns
    integer        :: nrows
!   integer,pointer:: n_row_node(:)=>null() ! # rows in each node
!   integer,pointer:: row_g2l(:)   =>null() ! Global to local row index
!   integer,pointer:: row_l2g(:)   =>null() ! Local to global row index
   integer,pointer:: n_col(:)     =>null() ! Nonzero cols of each row
   integer,pointer:: list_col(:)  =>null() ! Index of nonzero columns
   integer,pointer:: list_ptr(:)  =>null() ! First element of each row
 end type

 type sparseMatrix
   type(sparsity_t),pointer:: sparsity =>null() ! Nonzero-element indexes
   real(dp),pointer        :: val(:)   =>null() ! Nonzero-element values
 end type

CONTAINS

 subroutine copy_sparsity(a, b)
   type(sparsity_t), intent(in)  :: a
   type(sparsity_t), intent(out) :: b

   integer :: nrows, nnz

   nrows = a%nrows
   nnz = sum(a%n_col(:))

   b%nrows = nrows
   call re_alloc( b%n_col, 1,nrows)
   call re_alloc( b%list_ptr, 1,nrows)
   b%n_col = a%n_col
   b%list_ptr = a%list_ptr
   call re_alloc( b%list_col, 1,nnz)
   b%list_col = a%list_col

 end subroutine copy_sparsity



 subroutine setup_old_style_arrays(sm,nrows,num,listptr,list,val)
   type(sparseMatrix),     intent(in)  :: sm
   integer, intent(out) :: nrows
   integer, pointer     :: num(:), listptr(:)
   integer, pointer     :: list(:)
   real(dp), pointer     :: val(:)

   nrows    = sm%sparsity%nrows
   num     => sm%sparsity%n_col
   listptr => sm%sparsity%list_ptr
   list    => sm%sparsity%list_col
   val     => sm%val

 end subroutine setup_old_style_arrays

 subroutine build_sparsity(nrows,num,listptr,list,sparsity)
   type(sparsity_t),     intent(out)  :: sparsity
   integer, intent(in ) :: nrows
   integer, pointer     :: num(:), listptr(:)
   integer, pointer     :: list(:)
   real(dp), pointer     :: val(:)

   sparsity%nrows = nrows
   sparsity%n_col => num
   sparsity%list_ptr => listptr
   sparsity%list_col => list

 end subroutine build_sparsity

 subroutine sparse_alloc( matrix, sparsity )  ! Allocates a sparse matrix

   implicit none
   type(sparseMatrix),     intent(inout):: matrix
   type(sparsity_t),target,intent(in)   :: sparsity

   integer :: n_tot

   n_tot = sum(sparsity%n_col(:))
   call re_alloc( matrix%val, 1,n_tot, routine='sparse_alloc' )

   matrix%sparsity => sparsity

 end subroutine sparse_alloc

 subroutine copy_sparse_matrix ( a, b )
   type(sparseMatrix),     intent(in):: a
   type(sparseMatrix),     intent(out):: b

   integer nnz

   nnz = sum(a%sparsity%n_col(:))
   call re_alloc(b%val,1,nnz)
   b%val(:) = a%val(:)
   b%sparsity => a%sparsity

 end subroutine copy_sparse_matrix

 subroutine change_sparsity (m, new_sp, work)
   implicit none
   
   type(sparseMatrix),     intent(inout):: m
   type(sparsity_t),target,intent(in)   :: new_sp

   real(dp), pointer, optional :: work(:)

   integer :: size_old, size_new, i, in, ind, j, new_size
   integer :: maxval_j_old, maxval_j_new, max_col

   real(dp), pointer :: aux(:), tmp_val(:)
   type(sparsity_t), pointer :: old_sp

   old_sp => m%sparsity
   if (old_sp%nrows /= new_sp%nrows)  &
      call die("Incompatible nrows in change_sparsity")


   new_size = sum(new_sp%n_col(:))
   call re_alloc(tmp_val,1,new_size)
   tmp_val(:) = 0.0_dp

   maxval_j_old = maxval(old_sp%list_col(:))
   maxval_j_new = maxval(new_sp%list_col(:))
   max_col = max(maxval_j_old, maxval_j_new)

   if (present(work)) then
      call re_alloc(work,1,max_col)
      aux => work
   else
      nullify(aux)
      call re_alloc(aux,1,max_col)
   endif

   aux(1:max_col) = 0.0_dp
   do i=1,new_sp%nrows
      do in=1,old_sp%n_col(i)
         ind = old_sp%list_ptr(i) + in
         j = old_sp%list_col(ind)
         aux(j) = m%val(ind)
      enddo
      do in=1,new_sp%n_col(i)
         ind = new_sp%list_ptr(i) + in
         j = new_sp%list_col(ind)
         tmp_val(in) = aux(j)
      enddo
      aux(1:max_col) = 0.0_dp
   enddo

   call de_alloc(m%val)
   m%val => tmp_val

   ! The user is responsible for disposing of the old sparsity
   ! object when it is no longer needed

   !** AG DEALLOC m%sparsity (DEEP)
   m%sparsity => new_sp

   if (.not. present(work)) then
      call de_alloc(aux)
   endif

 end subroutine change_sparsity

END MODULE m_sparse_types
