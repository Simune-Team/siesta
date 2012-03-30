module class_SpMatrix

  use class_Sparsity
  use class_Array2D
  use class_OrbitalDistribution

  implicit none

  public :: SpMatrix
  public :: init, delete, assignment(=), refcount
  public :: val, spar, dist
  public :: nrows, nnzs, n_col, list_ptr, list_col
  public :: printSpMatrix
  public :: newSpMatrix

  private

  integer, parameter :: dp = selected_real_kind(10,100)

  ! This is the "meat" of the type
  !
  type SpMatrix_
    integer :: refCount = 0
    !----------------------
    character(len=256)   :: name = "null_SpMatrix"
    type(Sparsity)       :: sp
    type(Array2D)        :: a2d
    type(OrbitalDistribution)        :: dist
  end type SpMatrix_

  ! This is a wrapper type to be passed around
  type SpMatrix
    type(SpMatrix_), pointer :: data => null()
  end type SpMatrix

  interface assignment(=)
    module procedure assignSpMatrix
  end interface

  interface init
    module procedure initSpMatrix
  end interface

  interface delete
    module procedure deleteSpMatrix
  end interface

  interface newSpMatrix
    module procedure newSpMatrixFromArray2D
    module procedure newSpMatrixFromDims
  end interface

  interface refcount
     module procedure refcountSpMatrix
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

contains

   subroutine initSpMatrix(this)
     !........................................
     ! Constructor
     !........................................
     type (SpMatrix) :: this
     integer :: error

     call delete(this)
     allocate(this%data, stat=error)
     this%data%refCount = 1
     print *, "--> allocating SpMatrix "

  end subroutine initSpMatrix

  subroutine deleteSpMatrix(this)
    !............................
    ! Destructor
    !............................
    type (SpMatrix) :: this
    integer :: error

    if (.not. associated(this%data)) then
       ! print *, "--> deleting a non-associated SpMatrix"
       return
    endif
    print *, "--> attempting to delete SpMatrix:" // trim(this%data%name)
    this%data%refCount = this%data%refCount - 1
    if (this%data%refCount == 0) then
      ! Safe to delete the data now
      print *, "--> deallocating SpMatrix: " // trim(this%data%name) 
      call deleteSpMatrixData(this%data)
      deallocate(this%data, stat=error)
    endif

    this%data => null()

    CONTAINS
     subroutine deleteSpMatrixData(smdata)
      type(SpMatrix_) :: smdata

      call delete(smdata%sp)
      call delete(smdata%a2d)
      call delete(smdata%dist)
     end subroutine deleteSpMatrixData

  end subroutine deleteSpMatrix


  subroutine assignSpMatrix(this, other)
    !..................................................................
    ! Setting one list equal to another of the same type.
    ! No data copy, just increment reference and reference the same data
    !...................................................................
    type (SpMatrix), intent(inout) :: this
    type (SpMatrix), intent(in) :: other

    if (.not. associated(other%data)) then
     call die('Assignment to object that has not been initialized!')
    endif

    ! Delete to reset the pointers and decrement the ref count
    call delete(this)

    this%data => other%data
    this%data%refcount = this%data%refcount+1
    print *, "--> assigned SpMatrix: " // trim(this%data%name)

  end subroutine assignSpMatrix

  function refcountSpMatrix(this) result(count)
   type(SpMatrix), intent(in)  :: this
   integer :: count
   count = this%data%refCount
 end function refcountSpMatrix
    

  subroutine newSpMatrixFromArray2D(sp,a2d,dist,this,name)
     !........................................
     ! Constructor
     !........................................
     type (SpMatrix), intent(inout) :: this
     type(Sparsity), intent(in)   :: sp
     type(Array2D),  intent(in)   :: a2d
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
        this%data%name = "(from sp, dist, and a2d)"
     endif
     print *, "--> new SpMatrix: " // trim(this%data%name)

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
     call newArray2D(this%data%a2d,  &
                     nnzs(sp),dim2,"(new from SpMatrix)")

     if (present(name)) then
        this%data%name = trim(name)
     else
        this%data%name = "(from sp, dist, and dim2)"
     endif
     print *, "--> new SpMatrix: " // trim(this%data%name)

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
   call printSparsity(this%data%sp)
   call printArray2D(this%data%a2d)
   print "(a,i0,a)", "refcount: ",refcount(this),">"

 end subroutine printSpMatrix


 subroutine die(str)
  character(len=*), optional :: str
  if (present(str)) then
     print *, trim(str)
  endif
  stop
 end subroutine die

end module class_SpMatrix
