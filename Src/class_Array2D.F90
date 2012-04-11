module class_Array2D

  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)

  !
  type Array2D_
    integer :: refCount = 0
    !----------------------
    character(len=256)   :: name = "null Array2D"
    real(dp),   pointer  :: val(:,:) => null() ! Nonzero-element values
  end type Array2D_

  type Array2D
    type(Array2D_), pointer :: data => null()
  end type Array2D

  public  :: newArray2D, print_type, val

  interface assignment(=)
    module procedure Array2DfromNakedArray
  end interface

  interface val
    module procedure valArray2D
  end interface

  interface print_type
    module procedure printArray2D
  end interface

!========================
#define TYPE_NAME Array2D
#define __WHERE__ __FILE__
#include "basic_type.inc"
!========================

     subroutine delete_Data(a2d_data)
      type(Array2D_) :: a2d_data
      if (associated(a2d_data%val)) then
        deallocate( a2d_data%val)	
      endif
     end subroutine delete_Data


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
