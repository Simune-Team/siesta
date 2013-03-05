module class_Sparsity
  
  use alloc, only: re_alloc, de_alloc
  
  implicit none
  
  public :: newSparsity
  public :: print_type
  public :: nrows, nrows_g, nnzs, n_col, list_ptr, list_col

  character(len=*), parameter :: mod_name=__FILE__

  ! This is the "meat" of the type
  !
  type Sparsity_
    integer :: refCount = 0
    character(len=36)  :: id = "null_id"
    !----------------------
    character(len=256) :: name = "null_sparsity"
    integer            :: nrows = 0             ! Local number of rows
    integer            :: nrows_g = 0           ! Global number or rows
    integer            :: nnzs  = 0             ! Local number of non-zeros
    integer, pointer   :: n_col(:)     =>null() ! Nonzero cols of each row
    integer, pointer   :: list_col(:)  =>null() ! Index of nonzero columns
    integer, pointer   :: list_ptr(:)  =>null() ! First element of each row
    logical            :: initialized = .false.
  end type Sparsity_

  ! This is a wrapper type to be passed around
  type Sparsity
     type(Sparsity_), pointer :: data => null()
  end type Sparsity

  interface nrows
     module procedure nrowsSparsity
  end interface

  interface nrows_g
     module procedure nrows_gSparsity
  end interface

  interface nnzs
     module procedure nnzsSparsity
  end interface
  interface n_col
     module procedure n_colSparsity
     module procedure n_colSparsityI
  end interface
  interface list_ptr
     module procedure list_ptrSparsity
     module procedure list_ptrSparsityI
  end interface
  interface list_col
     module procedure list_colSparsity
     module procedure list_colSparsityI
  end interface

  interface print_type
     module procedure printSparsity
  end interface

!===========================
#define TYPE_NAME Sparsity
#include "basic_type.inc"
!===========================

  subroutine delete_Data(spdata)
    type(Sparsity_) :: spdata
    if (.not. spdata%initialized) RETURN
    call de_alloc( spdata%n_col,   &
         name="n_col " // trim(spdata%name),routine="Sparsity")	
    call de_alloc( spdata%list_ptr,   &
         name="list_ptr " // trim(spdata%name),routine="Sparsity")	
    call de_alloc( spdata%list_col,   &
         name="list_col " // trim(spdata%name),routine="Sparsity")	
  end subroutine delete_Data

!--------------------------------------------------------------------    
  subroutine newSparsity(sp,nrows,nrows_g,nnzs,num,listptr,list,name)

    type(Sparsity), intent(inout) :: sp

    integer, intent(in)           :: nrows, nrows_g, nnzs
    integer, intent(in)           :: num(:), listptr(:)
    integer, intent(in)           :: list(:)
    character(len=*), intent(in)  :: name

    integer :: stat

   ! We release the previous incarnation
   ! This means that we relinquish access to the previous
   ! memory location. It will be deallocated when nobody
   ! else is using it.

    call init(sp)
    
    sp%data%name = trim(name)
    
    call re_alloc(sp%data%n_col,1,nrows, &
         name="n_col " // trim(sp%data%name),routine="Sparsity")	
    call re_alloc(sp%data%list_ptr,1,nrows, &
         name="list_ptr " // trim(sp%data%name),routine="Sparsity")	

    sp%data%nrows = nrows
    sp%data%nrows_g = nrows_g
    sp%data%nnzs  = nnzs
    sp%data%n_col(1:nrows) = num(1:nrows)
    sp%data%list_ptr(1:nrows) = listptr(1:nrows)

    if (nnzs /= sum(num(1:nrows))) then
       call die("nnzs mismatch in new_sparsity")
    endif

    call re_alloc(sp%data%list_col,1,nnzs, &
         name="list_col " // trim(sp%data%name),routine="Sparsity")	
    sp%data%list_col(1:nnzs) = list(1:nnzs)

    sp%data%initialized = .true.   

    call tag_new_object(sp)

  end subroutine newSparsity
  
  function nrowsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nrows
  end function nrowsSparsity

  function nrows_gSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nrows_g
  end function nrows_gSparsity

  function nnzsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nnzs
  end function nnzsSparsity

  function n_colSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%n_col
  end function n_colSparsity
  function n_colSparsityI(this,i) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: i
    integer                    :: p
    p = this%data%n_col(i)
  end function n_colSparsityI


  function list_ptrSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%list_ptr
  end function list_ptrSparsity
  function list_ptrSparsityI(this,i) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: i
    integer                    :: p
    p = this%data%list_ptr(i)
  end function list_ptrSparsityI

  function list_colSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%list_col
  end function list_colSparsity
  function list_colSparsityI(this,i) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: i
    integer                    :: p
    p = this%data%list_col(i)
  end function list_colSparsityI

  subroutine printSparsity(sp)
    type(Sparsity), intent(in) :: sp

    if (.not. associated(sp%data)) then
       print "(a)", "Sparsity Not Associated"
       RETURN
    endif

    print "(a,i0,a,i0,a,i0,a,i0,a)", &
                "  <sparsity:"//trim(sp%data%name)//" nrows_g=", &
                sp%data%nrows_g, " nrows=", sp%data%nrows, &
                " nnzs=",sp%data%nnzs,", refcount: ",  &
                refcount(sp), ">"
  end subroutine printSparsity

end module class_Sparsity
