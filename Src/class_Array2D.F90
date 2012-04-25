module class_Array2D

  use alloc, only: re_alloc, de_alloc

  implicit none

  character(len=*), parameter :: mod_name=__FILE__

  integer, parameter :: dp = selected_real_kind(10,100)

  !
  type Array2D_
    integer :: refCount = 0
    character(len=36)   :: id = "null_id"
    !----------------------
    character(len=256)   :: name = "null Array2D"
    real(dp),   pointer  :: val(:,:) => null() ! Nonzero-element values
  end type Array2D_

  type Array2D
    type(Array2D_), pointer :: data => null()
  end type Array2D

  public  :: newArray2D, print_type, val

  interface newArray2D
    module procedure newArray2DfromDimensions
    module procedure newArray2DfromNakedArray
  end interface

  interface val
    module procedure valArray2D
  end interface

  interface print_type
    module procedure printArray2D
  end interface

!========================
#define TYPE_NAME Array2D
#include "basic_type.inc"
!========================

     subroutine delete_Data(a2d_data)
      type(Array2D_) :: a2d_data
      if (associated(a2d_data%val)) then
        call de_alloc( a2d_data%val, &
             name="val "//trim(a2d_data%name),routine="Array2D")	
      endif
     end subroutine delete_Data


  subroutine newArray2DFromDimensions(this,n,m,name)
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

   if (present(name)) then
      this%data%name = trim(name)
   else
      this%data%name = "(Array2D from n,m)"
   endif

   call re_alloc(this%data%val,1,n,1,m, &
        name="val "//trim(this%data%name),routine="Array2D")	
   this%data%val(:,:) = 0.0_dp

   call tag_new_object(this)

 end subroutine newArray2DFromDimensions

  subroutine newArray2DfromNakedArray(this, val, name)
    !..................................................................
    !...................................................................
    type (Array2D), intent(inout) :: this
    real(dp), intent(in) :: val(:,:)
    character(len=*), intent(in), optional  :: name

    integer :: n, m

    call init(this)

    n = size(val,dim=1)
    m = size(val,dim=2)

    if (present(name)) then
       this%data%name = trim(name)
    else
       this%data%name = "(Array2D from naked array)"
    endif

    call re_alloc(this%data%val,1,n,1,m, &
        name="val "//trim(this%data%name),routine="Array2D")	
    this%data%val(:,:) = val(:,:)

   call tag_new_object(this)

 end subroutine newArray2DfromNakedArray


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

end module class_Array2D
