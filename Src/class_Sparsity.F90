module class_Sparsity
  
  use alloc, only: re_alloc, de_alloc
  
  implicit none
  
  public :: newSparsity
  public :: print_type
  public :: attach
  public :: nrows, nrows_g
  public :: ncols, ncols_g
  public :: n_row, n_col
  public :: nnzs, list_ptr, list_col

  character(len=*), parameter :: mod_name= __FILE__

  ! This is the "meat" of the type
  !
  type Sparsity_
    integer :: refCount = 0
    character(len=36)  :: id = "null_id"
    !----------------------
    character(len=256) :: name = "null_sparsity"
    integer            :: nrows = 0             ! Local number of rows
    integer            :: nrows_g = 0           ! Global number or rows
    integer            :: ncols = 0             ! Local number of columns
    integer            :: ncols_g = 0           ! Global number or columns
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

  interface ncols
     module procedure ncolsSparsity
  end interface

  interface ncols_g
     module procedure ncols_gSparsity
  end interface

! ******************
! Specific routines to retrieve direct information
! about the sparsity entries.

  interface attach
     module procedure attachSparsity
  end interface

  interface n_col
     module procedure n_colSparsity
     module procedure n_colSparsityI
  end interface

  interface n_row
     ! This interface ensures that one can retrieve
     ! all information from the sparse class.
     module procedure n_rowSparsityI
  end interface

  interface nnzs
     module procedure nnzsSparsity
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
! For easy compatability we have added ncols and ncols_g as
! optional arguments.
! In case, not specified, values of nrows_g are used. (block-cyclic)
  subroutine newSparsity(sp,nrows,nrows_g,nnzs,num,listptr,list,name, &
       ncols,ncols_g)

    type(Sparsity), intent(inout) :: sp

    integer, intent(in)           :: nrows, nrows_g, nnzs
    integer, intent(in)           :: num(:), listptr(:)
    integer, intent(in)           :: list(:)
    character(len=*), intent(in)  :: name
    integer, intent(in), optional :: ncols, ncols_g

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
    if ( present(ncols_g) ) then
       sp%data%ncols_g = ncols_g
    else
       ! We default it to a Block-cyclic distribution, hence
       ! the number of columns is the same as the number of 
       ! global rows (we cannot guess super-cells, that would indeed be amazing)
       sp%data%ncols_g = nrows_g
    end if
    if ( present(ncols) ) then
       sp%data%ncols = ncols
    else
       ! Again, block cyclic has the maximum number of columns
       ! in each block
       sp%data%ncols = sp%data%ncols_g
    end if
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
  
  pure function nrowsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nrows
  end function nrowsSparsity

  elemental function n_rowSparsityI(this,col) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: col
    integer                    :: p
    p = count(this%data%list_col == col)
  end function n_rowSparsityI

  pure function nrows_gSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nrows_g
  end function nrows_gSparsity

  pure function ncolsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%ncols
  end function ncolsSparsity

  pure function ncols_gSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%ncols_g
  end function ncols_gSparsity


  pure function nnzsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nnzs
  end function nnzsSparsity

  function n_colSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%n_col
  end function n_colSparsity
  elemental function n_colSparsityI(this,row) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: row
    integer                    :: p
    p = this%data%n_col(row)
  end function n_colSparsityI

  function list_ptrSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%list_ptr
  end function list_ptrSparsity
  elemental function list_ptrSparsityI(this,i) result (p)
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
  elemental function list_colSparsityI(this,i) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: i
    integer                    :: p
    p = this%data%list_col(i)
  end function list_colSparsityI

  ! Generic routine for retrieval of information in one line...
  subroutine attachSparsity(this,D,n_col,list_col,list_ptr, &
       nrows,nrows_g,ncols,ncols_g, &
       nnzs)
    type(Sparsity), intent(inout) :: this
    logical, intent(in) , optional :: D ! DUMMY, force named arguments
    ! Optional arguments to retrieve all information with named arguments
    integer, pointer    , optional :: n_col(:), list_col(:), list_ptr(:)
    integer, intent(out), optional :: nrows, nrows_g
    integer, intent(out), optional :: ncols, ncols_g
    integer, intent(out), optional :: nnzs
    if ( present(D) ) call die('PROGRAMMING ERROR, named args please')

    ! the arrays
    if ( present(n_col) ) n_col => n_colSparsity(this)
    if ( present(list_col) ) list_col => list_colSparsity(this)
    if ( present(list_ptr) ) list_ptr => list_ptrSparsity(this)
    ! the integers
    if ( present(nrows) ) nrows = nrowsSparsity(this)
    if ( present(nrows_g) ) nrows_g = nrows_gSparsity(this)
    if ( present(ncols) ) ncols = ncolsSparsity(this)
    if ( present(ncols_g) ) ncols_g = ncols_gSparsity(this)
    if ( present(nnzs) ) nnzs = nnzsSparsity(this)

  end subroutine attachSparsity

  subroutine printSparsity(sp)
    type(Sparsity), intent(in) :: sp

    if (.not. initialized(sp) ) then
       print "(a)", "Sparsity Not Associated"
       RETURN
    endif

    print "(2(a,i0),a,f0.4,2(a,i0),a)", &
                "  <sparsity:"//trim(sp%data%name)//" nrows_g=", &
                sp%data%nrows_g, " nrows=", sp%data%nrows, &
                " sparsity=",real(sp%data%nnzs)/ &
                (sp%data%nrows_g*sp%data%ncols_g), &
                " nnzs=",sp%data%nnzs,", refcount: ",  &
                refcount(sp), ">"
  end subroutine printSparsity

end module class_Sparsity
