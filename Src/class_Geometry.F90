module class_Geometry

  implicit none

  character(len=*), parameter :: mod_name=__FILE__

  integer, parameter :: dp = selected_real_kind(10,100)

  !
  type Geometry_
    integer            :: refCount = 0
    character(len=36)  :: id = "null_id"
    !----------------------
    character(len=256)   :: name = "null Geometry"
    integer              :: na                ! Number of atoms
    real(dp)             :: cell(3,3)         ! Cell vectors by columns
    real(dp),   pointer  :: xa(:,:) => null() ! Nonzero-element values
    integer,    pointer  :: isa(:) => null()  ! Species list
  end type Geometry_

  type Geometry
    type(Geometry_), pointer :: data => null()
  end type Geometry

  public  :: newGeometry, print_type, coords

  interface coords
    module procedure coordsGeometry
  end interface

  interface print_type
    module procedure printGeometry
  end interface

!========================
#define TYPE_NAME Geometry
#include "basic_type.inc"
!========================

     subroutine delete_Data(a2d_data)
      type(Geometry_) :: a2d_data
      if (associated(a2d_data%xa)) then
        deallocate( a2d_data%xa)
      endif
      if (associated(a2d_data%isa)) then
        deallocate( a2d_data%isa)
      endif
     end subroutine delete_Data


  subroutine newGeometry(this,na,cell,xa,isa,name)

   type(Geometry), intent(inout)  :: this
   integer, intent(in)            :: na
   real(dp), intent(in)           :: cell(3,3)
   real(dp), intent(in)           :: xa(3,na)
   integer, intent(in)            :: isa(na)
   character(len=*), intent(in), optional  :: name

   ! We release the previous incarnation
   ! This means that we relinquish access to the previous
   ! memory location. It will be deallocated when nobody
   ! else is using it.

   call init(this)

   allocate(this%data%xa(3,na))
   allocate(this%data%isa(na))

   this%data%na = na
   this%data%cell(:,:) = cell(:,:)
   this%data%xa(:,:) = xa(:,:)
   this%data%isa(:) = isa(:)

   if (present(name)) then
      this%data%name = trim(name)
   else
      this%data%name = "(Geometry)"
   endif
   call tag_new_object(this)

  end subroutine newGeometry

  function coordsGeometry(this) result(p)
   type(Geometry), intent(in)  :: this
   real(dp), pointer           :: p(:,:)

   nullify(p)
   p => this%data%xa
 end function coordsGeometry

 subroutine printGeometry(this)
   type(Geometry), intent(in)  :: this

   if (.not. associated(this%data)) then
      print "(a)", "Geometry Not Associated"
      RETURN
   endif

   print "(a,i0,a,i0,a)", "  <Geometry:" // trim(this%data%name) // &
                               " na=",  this%data%na,  &
                               ", refcount: ", refcount(this),">"
 end subroutine printGeometry


 subroutine die(str)
  character(len=*), optional :: str
  if (present(str)) then
     print *, trim(str)
  endif
  stop
 end subroutine die

end module class_Geometry
