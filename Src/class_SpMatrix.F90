module class_SpMatrix

  use class_Sparsity
  use class_dArray2D
  use class_OrbitalDistribution

  implicit none

  character(len=*), parameter :: mod_name=__FILE__

  public :: val, spar, dist
  public :: nrows, nrows_g, nnzs, n_col, list_ptr, list_col
  public :: print_type
  public :: newSpMatrix

  integer, parameter :: dp = selected_real_kind(10,100)

  type SpMatrix_
    integer            :: refCount = 0
    character(len=36)  :: id = "null_id"
    !----------------------
    character(len=256)   :: name = "null_SpMatrix"
    type(Sparsity)       :: sp
    type(dArray2D)       :: a2d
    type(OrbitalDistribution)        :: dist
  end type SpMatrix_

  type SpMatrix
    type(SpMatrix_), pointer :: data => null()
  end type SpMatrix

  interface newSpMatrix
    module procedure newSpMatrixFromArray2D
    module procedure newSpMatrixFromDims
  end interface

  interface val
    module procedure valSpMatrix
  end interface

  interface spar
    module procedure sparSpMatrix
  end interface

  interface dist
    module procedure distSpMatrix
  end interface

  interface nrows
     module procedure nrowsSpMatrix
  end interface

  interface nrows_g
     module procedure nrows_gSpMatrix
  end interface

  interface nnzs
     module procedure nnzsSpMatrix
  end interface

  interface n_col
     module procedure n_colSpMatrix
  end interface

  interface list_ptr
     module procedure list_ptrSpMatrix
  end interface

  interface list_col
     module procedure list_colSpMatrix
  end interface

interface print_type
   module procedure printSpMatrix
end interface

!==========================
#define TYPE_NAME SpMatrix
#include "basic_type.inc"
!==========================

     subroutine delete_Data(smdata)
      type(SpMatrix_) :: smdata

      call delete(smdata%sp)
      call delete(smdata%a2d)
      call delete(smdata%dist)
    end subroutine delete_Data

  subroutine newSpMatrixFromArray2D(sp,a2d,dist,this,name)
     !........................................
     ! Constructor
     !........................................
     type (SpMatrix), intent(inout) :: this
     type(Sparsity), intent(in)   :: sp
     type(dArray2D),  intent(in)  :: a2d
     type(OrbitalDistribution),  intent(in)   :: dist
     character(len=*), intent(in), optional :: name

     integer :: error

     call init(this)

     this%data%sp = sp
     this%data%a2d = a2d
     this%data%dist = dist

     if (present(name)) then
        this%data%name = trim(name)
     else
        this%data%name = "(SpMatrix from sp, dist, and a2d)"
     endif
     call tag_new_object(this)

   end subroutine newSpMatrixFromArray2D

  subroutine newSpMatrixFromDims(sp,dim2,dist,this,name)
     !........................................
     ! Constructor
     !........................................
     type (SpMatrix), intent(inout) :: this
     type(Sparsity), intent(in)   :: sp
     type(OrbitalDistribution), intent(in)   :: dist
     integer,  intent(in)         :: dim2
     character(len=*), intent(in), optional :: name

     integer :: error

     call init(this)
     this%data%sp = sp
     this%data%dist = dist
     call newdArray2D(this%data%a2d,  &
                     nnzs(sp),dim2,"(new from SpMatrix)")

     if (present(name)) then
        this%data%name = trim(name)
     else
        this%data%name = "(SpMatrix from sp, dim2, and dist)"
     endif
     call tag_new_object(this)

   end subroutine newSpMatrixFromDims

!--------------------------------------------------
  function valSpMatrix(this) result(p)
   type(SpMatrix), intent(in)  :: this
   real(dp), pointer          :: p(:,:)

   nullify(p)
   p => val(this%data%a2d)
 end function valSpMatrix

  function sparSpMatrix(this) result(p)
   type(SpMatrix), intent(in)  :: this
   type(Sparsity), pointer     :: p

   nullify(p)
   p => this%data%sp
 end function sparSpMatrix

  function distSpMatrix(this) result(p)
   type(SpMatrix), intent(in)  :: this
   type(OrbitalDistribution), pointer     :: p

   nullify(p)
   p => this%data%dist
 end function distSpMatrix

!--------------------------------------------------
 function nrowsSpmatrix(this) result (n)
   type(SpMatrix), intent(in)  :: this
   integer                     :: n
   n = nrows(this%data%sp)
 end function nrowsSpmatrix

 function nrows_gSpmatrix(this) result (n)
   type(SpMatrix), intent(in)  :: this
   integer                     :: n
   n = nrows_g(this%data%sp)
 end function nrows_gSpmatrix

 function nnzsSpMatrix(this) result (n)
   type(SpMatrix), intent(in)  :: this
   integer                     :: n
   n = nnzs(this%data%sp)
 end function nnzsSpMatrix

 function n_colSpMatrix(this) result (p)
   type(SpMatrix), intent(in)  :: this
   integer, pointer            :: p(:)
   p => n_col(this%data%sp)
 end function n_colSpMatrix

 function list_ptrSpMatrix(this) result (p)
   type(SpMatrix), intent(in)  :: this
   integer, pointer            :: p(:)
   p => list_ptr(this%data%sp)
 end function list_ptrSpMatrix

 function list_colSpMatrix(this) result (p)
   type(SpMatrix), intent(in)  :: this
   integer, pointer            :: p(:)
   p => list_col(this%data%sp)
 end function list_colSpMatrix
   
 subroutine printSpMatrix(this)
   type(SpMatrix), intent(in)  :: this

   if (.not. associated(this%data)) then
      print "(a)", "SpMatrix Not Associated"
      RETURN
   endif

   print "(a)", "<spMatrix:" // trim(this%data%name)
   call print_type(this%data%sp)
   call print_type(this%data%a2d)
   print "(a,i0,a)", "refcount: ",refcount(this),">"

 end subroutine printSpMatrix

end module class_SpMatrix
