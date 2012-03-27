module class_Sparsity

  implicit none

  public :: Sparsity
  public :: init, delete, assignment(=), refcount
  private :: assignSparsity, initSparsity, deleteSparsity, refcountSparsity
  public :: newSparsity
  public :: printSparsity
  public :: nrows, nnzs, n_col, list_ptr, list_col
  private :: nrowsSparsity, nnzsSparsity, n_colSparsity
  private :: list_ptrSparsity, list_colSparsity
  private

  ! This is the "meat" of the type
  !
  type Sparsity_
    integer :: refCount = 0
    !----------------------
    character(len=256) :: name = "null_sparsity"
!   integer        :: n_row_g=0  ! Global number of rows
!   integer        :: n_col_g=0  ! Global number of columns
    integer        :: nrows = 0
    integer        :: ncols = 0
    integer        :: nnzs  = 0
!   integer,pointer:: n_row_node(:)=>null() ! # rows in each node
!   integer,pointer:: row_g2l(:)   =>null() ! Global to local row index
!   integer,pointer:: row_l2g(:)   =>null() ! Local to global row index
   integer,pointer:: n_col(:)     =>null() ! Nonzero cols of each row
   integer,pointer:: list_col(:)  =>null() ! Index of nonzero columns
   integer,pointer:: list_ptr(:)  =>null() ! First element of each row
   logical        :: initialized = .false.
  end type Sparsity_

  ! This is a wrapper type to be passed around
  type Sparsity
    type(Sparsity_), pointer :: data => null()
  end type Sparsity

  interface assignment(=)
    module procedure assignSparsity
  end interface

  interface init
    module procedure initSparsity
  end interface

  interface delete
    module procedure deleteSparsity
  end interface

  interface refcount
    module procedure refcountSparsity
  end interface

  interface nrows
     module procedure nrowsSparsity
  end interface
  interface nnzs
     module procedure nnzsSparsity
  end interface
  interface n_col
     module procedure n_colSparsity
  end interface
  interface list_ptr
     module procedure list_ptrSparsity
  end interface
  interface list_col
     module procedure list_colSparsity
  end interface

contains

   subroutine initSparsity(this)
     !........................................
     ! Constructor
     !........................................
     type (Sparsity) :: this
     integer :: error

     call delete(this)
     allocate(this%data, stat=error)
     this%data%refCount = 1
     print *, "--> allocated sparsity "

  end subroutine initSparsity

  subroutine deleteSparsity(this)
    !............................
    ! Destructor
    !............................
    type (Sparsity) :: this
    integer :: error

    if (.not. associated(this%data)) return
    print *, "... attempting to deallocate sparsity: " //  &        
                           trim(this%data%name) 

    this%data%refCount = this%data%refCount - 1
    if (this%data%refCount == 0) then
      ! Safe to delete the data now
      call deleteSparsityData(this%data)
      print *, "--> deallocated sparsity: " // trim(this%data%name) 
      deallocate(this%data, stat=error)
    endif

    this%data => null()

    CONTAINS
     subroutine deleteSparsityData(spdata)
      type(Sparsity_) :: spdata
      if (.not. spdata%initialized) RETURN
      deallocate( spdata%n_col)
      deallocate( spdata%list_ptr)
      deallocate( spdata%list_col)
     end subroutine deleteSparsityData

  end subroutine deleteSparsity


  subroutine assignSparsity(this, other)
    !..................................................................
    ! Setting one list equal to another of the same type.
    ! No data copy, just increment reference and reference the same data
    !...................................................................
    type (Sparsity), intent(inout) :: this
    type (Sparsity), intent(in) :: other

    if (.not. associated(other%data)) then
     call die('Assignment to object that has not been initialized!')
    endif

    ! Delete to reset the pointers and decrement the ref count
    call delete(this)

    this%data => other%data
    this%data%refcount = this%data%refcount+1

  end subroutine assignSparsity

  function refcountSparsity(this) result(count)
   type(Sparsity), intent(in)  :: this
   integer :: count
   count = this%data%refCount
  end function refcountSparsity
    
 subroutine newSparsity(sp,nrows,ncols,nnzs,num,listptr,list,name)

   type(Sparsity), intent(inout)  :: sp

   integer, intent(in)  :: nrows, ncols, nnzs
   integer, intent(in)  :: num(:), listptr(:)
   integer, intent(in)  :: list(:)
   character(len=*), intent(in) :: name

   integer :: stat

   ! We release the previous incarnation
   ! This means that we relinquish access to the previous
   ! memory location. It will be deallocated when nobody
   ! else is using it.

   call init(sp)

   allocate(sp%data%n_col(1:nrows))
   allocate(sp%data%list_ptr(1:nrows))

   sp%data%nrows = nrows
   sp%data%ncols = ncols
   sp%data%nnzs  = nnzs
   sp%data%n_col(1:nrows) = num(1:nrows)
   sp%data%list_ptr(1:nrows) = listptr(1:nrows)

   if (nnzs /= sum(num(1:nrows))) then
      call die("nnzs mismatch in new_sparsity")
   endif

   allocate(sp%data%list_col(1:nnzs))
   sp%data%list_col(1:nnzs) = list(1:nnzs)

   sp%data%initialized = .true.   
   sp%data%name = trim(name)

   print *, "--> allocated sparsity data: " // trim(sp%data%name)
   
 end subroutine newSparsity

 function nrowsSparsity(this) result (n)
   type(Sparsity), intent(in)  :: this
   integer                     :: n
   n = this%data%nrows
 end function nrowsSparsity

 function nnzsSparsity(this) result (n)
   type(Sparsity), intent(in)  :: this
   integer                     :: n
   n = this%data%nnzs
 end function nnzsSparsity

 function n_colSparsity(this) result (p)
   type(Sparsity), intent(in)  :: this
   integer, pointer            :: p(:)
   p => this%data%n_col
 end function n_colSparsity

 function list_ptrSparsity(this) result (p)
   type(Sparsity), intent(in)  :: this
   integer, pointer            :: p(:)
   p => this%data%list_ptr
 end function list_ptrSparsity

 function list_colSparsity(this) result (p)
   type(Sparsity), intent(in)  :: this
   integer, pointer            :: p(:)
   p => this%data%list_col
 end function list_colSparsity
   
 subroutine die(str)
  character(len=*), optional :: str
  if (present(str)) then
     print *, trim(str)
  endif
  stop
 end subroutine die

 subroutine printSparsity(sp)
   type(Sparsity), intent(in)  :: sp

   if (.not. associated(sp%data)) then
      print "(a)", "Sparsity Not Associated"
      RETURN
   endif

   print "(a,i0,a,i0,a,i0,a)",     &
                          "  <sparsity:" // trim(sp%data%name) // " nrows=", &
                          sp%data%nrows, &
                          " nnzs=",sp%data%nnzs,", refcount: ",  &
                          refcount(sp), ">"
 end subroutine printSparsity

end module class_Sparsity
