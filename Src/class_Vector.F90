! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module class_Vector

  use alloc, only: re_alloc, de_alloc

  implicit none

  character(len=*), parameter :: mod_name= "Vector"

  integer, parameter :: dp = selected_real_kind(10,100)

  !
  type Vector_
    integer :: refCount = 0
    character(len=36)   :: id = "null_id"
    !----------------------
    character(len=256)   :: name = "null Vector"
    real(dp),   pointer  :: val(:) => null() ! Nonzero-element values
  end type Vector_

  type Vector
    type(Vector_), pointer :: data => null()
  end type Vector

  public  :: newVector, print_type, val

  interface newVector
    module procedure newVectorfromDimension
    module procedure newVectorfromNakedArray
  end interface

  interface val
    module procedure valVector
  end interface

  interface print_type
    module procedure printVector
  end interface

!========================
#define TYPE_NAME Vector
#include "basic_type.inc"
!========================

     subroutine delete_Data(vec_data)
      type(Vector_) :: vec_data
      if (associated(vec_data%val)) then
        call de_alloc( vec_data%val, &
             name="val "//trim(vec_data%name),routine="Vector")	
      endif
     end subroutine delete_Data


  subroutine newVectorFromDimension(this,n,name)
  ! This could be implemented also as an assignment 
  ! (see below)

   type(Vector), intent(inout)  :: this
   integer, intent(in)           :: n
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
      this%data%name = "(Vector from n)"
   endif

   call re_alloc(this%data%val,1,n, &
        name="val "//trim(this%data%name),routine="Vector")	
   this%data%val(:) = 0.0_dp

   call tag_new_object(this)

 end subroutine newVectorFromDimension

  subroutine newVectorfromNakedArray(this, val, name)
    !..................................................................
    !...................................................................
    type (Vector), intent(inout) :: this
    real(dp), intent(in) :: val(:)
    character(len=*), intent(in), optional  :: name

    integer :: n

    call init(this)

    n = size(val,dim=1)

    if (present(name)) then
       this%data%name = trim(name)
    else
       this%data%name = "(Vector from naked array)"
    endif

    call re_alloc(this%data%val,1,n, &
        name="val "//trim(this%data%name),routine="Vector")	
    this%data%val(:) = val(:)

   call tag_new_object(this)

 end subroutine newVectorfromNakedArray


  function valVector(this) result(p)
   type(Vector), intent(in)  :: this
   real(dp), pointer          :: p(:)

   nullify(p)
   p => this%data%val
 end function valVector

 subroutine printVector(this)
   type(Vector), intent(in)  :: this

    integer :: n, m

   if (.not. associated(this%data)) then
      print "(a)", "Vector Not Associated"
      RETURN
   endif

    n = size(this%data%val,dim=1)

   print "(a,i0,a,i0,a,i0,a)", "  <array:" // trim(this%data%name) // &
                               " n=",  n,   &
                               ", refcount: ", refcount(this),">"
 end subroutine printVector

end module class_Vector
