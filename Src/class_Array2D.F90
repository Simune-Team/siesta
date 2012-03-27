module class_Array2D

  implicit none

  public  :: Array2D
  public  :: init, delete, assignment(=), refcount
  private :: assignArray2D, initArray2D, deleteArray2D, refcountArray2D
  private :: Array2DfromNakedArray
  public  :: newArray2D
  public  :: printArray2D
  public  :: val
  private :: valArray2D
  private

  integer, parameter :: dp = selected_real_kind(10,100)

  ! This is the "meat" of the type
  !
  type Array2D_
    integer :: refCount = 0
    !----------------------
    character(len=256)   :: name = "null_array"
    real(dp),   pointer  :: val(:,:) => null() ! Nonzero-element values
  end type Array2D_

  ! This is a wrapper type to be passed around
  type Array2D
    type(Array2D_), pointer :: data => null()
  end type Array2D

  interface assignment(=)
    module procedure assignArray2D
    module procedure Array2DfromNakedArray
  end interface

  interface init
    module procedure initArray2D
  end interface

  interface delete
    module procedure deleteArray2D
  end interface

  interface refcount
    module procedure refcountArray2D
  end interface

  interface val
    module procedure valArray2D
  end interface

contains

   subroutine initArray2D(this)
     !........................................
     ! Constructor
     !........................................
     type (Array2D) :: this
     integer :: error

     call delete(this)
     allocate(this%data, stat=error)
     print *, "--> allocated 2DArray: " // trim(this%data%name) 
     this%data%refCount = 1

  end subroutine initArray2D

  subroutine deleteArray2D(this)
    !............................
    ! Destructor
    !............................
    type (Array2D) :: this
    integer :: error

    if (.not. associated(this%data)) return

    this%data%refCount = this%data%refCount - 1
    if (this%data%refCount == 0) then
      ! Safe to delete the data now
      print *, "--> deallocated 2DArray: " // trim(this%data%name) 
      call deleteArray2DData(this%data)
      deallocate(this%data, stat=error)
    endif

    this%data => null()

    CONTAINS
     subroutine deleteArray2DData(a2d_data)
      type(Array2D_) :: a2d_data
      if (associated(a2d_data%val)) then
        deallocate( a2d_data%val)	
      endif
     end subroutine deleteArray2DData

  end subroutine deleteArray2D


  subroutine assignArray2D(this, other)
    !..................................................................
    ! Setting one list equal to another of the same type.
    ! No data copy, just increment reference and reference the same data
    !...................................................................
    type (Array2D), intent(inout) :: this
    type (Array2D), intent(in) :: other

    if (.not. associated(other%data)) then
     call die('Assignment to object that has not been initialized!')
    endif

    ! Delete to reset the pointers and decrement the ref count
    call delete(this)

    this%data => other%data
    this%data%refcount = this%data%refcount+1

  end subroutine assignArray2D

  subroutine newArray2D(this,n,m,name)
  ! This could be implemented also as an assignment 
  ! (see below)

   type(Array2D), intent(inout)  :: this
   integer, intent(in)           :: n, m
   character(len=*), intent(in), optional  :: name

   integer :: stat

   ! We release the previous incarnation
   ! This means that we relinquish access to the previous
   ! memory location. It will be deallocated when nobody
   ! else is using it.

   call init(this)

   allocate(this%data%val(1:n,1:m))
   this%data%val(:,:) = 0.0_dp
   if (present(name)) then
      this%data%name = trim(name)
   else
      this%data%name = "(built from n,m)"
   endif

  end subroutine newArray2D

  subroutine Array2DfromNakedArray(this, val)
    !..................................................................
    !...................................................................
    type (Array2D), intent(inout) :: this
    real(dp), intent(in) :: val(:,:)

    integer :: n, m

    call init(this)

   n = size(val,dim=1)
   m = size(val,dim=2)

   allocate(this%data%val(1:n,1:m))
   this%data%val(:,:) = val(:,:)
   this%data%name = "(copied from naked array)"
  end subroutine Array2DfromNakedArray


  function refcountArray2D(this) result(count)
   type(Array2D), intent(in)  :: this
   integer  :: count
   count = this%data%refCount
  end function refcountArray2D

  function valArray2D(this) result(p)
   type(Array2D), intent(in)  :: this
   real(dp), pointer          :: p(:,:)

   nullify(p)
   p => this%data%val
 end function valArray2D

 subroutine printArray2D(this)
   type(Array2D), intent(in)  :: this

    integer :: n, m

   if (.not. associated(this%data)) then
      print "(a)", "Array2D Not Associated"
      RETURN
   endif

    n = size(this%data%val,dim=1)
    m = size(this%data%val,dim=2)

   print "(a,i0,a,i0,a,i0,a)", "  <array2D:" // trim(this%data%name) // &
                               " n=",  n," m=",m,   &
                               ", refcount: ", refcount(this),">"
 end subroutine printArray2D


 subroutine die(str)
  character(len=*), optional :: str
  if (present(str)) then
     print *, trim(str)
  endif
  stop
 end subroutine die

end module class_Array2D
