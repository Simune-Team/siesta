!--------------------------------------------------------------
!
! Basic object functionality
!
!   This file has to be included after the global declarations
!   for extra functionality, just where a "contains" statement
!   would go
!
!   Parameter:  "TYPE_NAME" has to be set to the name of the type
!--------------------------------------------------------------

  PRIVATE

  public  :: TYPE_NAME
  public  :: init, delete, assignment(=), refcount, id
  public  :: name
  public  :: same ! same %data address (i.e. if: this1 = this2)
  public  :: initialized ! is allocated

  interface assignment(=)
    module procedure assign_
  end interface

  interface init
    module procedure init_
  end interface

  interface delete
    module procedure delete_
  end interface

  interface refcount
    module procedure refcount_
  end interface

  interface id
    module procedure id_
  end interface

  interface name
    module procedure name_
  end interface

  interface initialized
    module procedure initialized_
  end interface

  interface same
    module procedure same_
  end interface


  ! Stand-alone routine which must be provided
  interface 
      subroutine die(str)
      character(len=*), intent(in)  :: str
      end subroutine die
  end interface

CONTAINS

! -- Main structural features
!
   subroutine init_(this)

     type (TYPE_NAME), intent(inout) :: this

     integer :: error

     call delete(this)
     allocate(this%data, stat=error)
     if (error /= 0) then
       call die("Allocation error!")
     endif
     this%data%refCount = 1

  end subroutine init_

  subroutine delete_(this)

    type (TYPE_NAME), intent(inout) :: this

    integer :: error
    logical, external  :: print_debug_object_info

    if (.not. initialized(this) ) return

    this%data%refCount = this%data%refCount - 1
    if (this%data%refCount == 0) then
      ! Safe to delete the data now
      call delete_Data(this%data)
      if (print_debug_object_info()) then
        print *, "--> deallocated " // id(this) // " " // trim(this%data%name) 
      endif
      deallocate(this%data, stat=error)
      if (error /= 0) then
         call die("Deallocation error!")
      endif
    endif

    ! This is important to use the correct initialized functions
    this%data => null()

  end subroutine delete_


  subroutine assign_(this, other)
    !..................................................................
    ! Setting one list equal to another of the same type.
    ! No data copy, just increment reference and reference the same data
    !...................................................................
    type (TYPE_NAME), intent(inout) :: this
    type (TYPE_NAME), intent(in) :: other

    if (.not. initialized(other) ) then
     call die("Assignment of non-initialized object in " // trim(mod_name))
    endif

    ! Delete to reset the pointers and decrement the ref count
    call delete(this)

    this%data => other%data
    this%data%refcount = this%data%refcount+1

  end subroutine assign_

  pure function initialized_(this) result(init)
    type(TYPE_NAME), intent(in) :: this
    logical :: init
    ! If it is not associated, it can not contain any data
    init = associated(this%data)
  end function initialized_

  pure function same_(this1,this2) result(same)
    type(TYPE_NAME), intent(in) :: this1, this2
    logical :: same
    ! If they are not both initialized they can not be the same
    same = initialized(this1) .and. initialized(this2)
    if ( .not. same ) return
    same = associated(this1%data, target=this2%data)
  end function same_

  function refcount_(this) result(count)
   type(TYPE_NAME), intent(in)  :: this
   integer  :: count
   count = this%data%refCount
  end function refcount_

  function id_(this) result(str)
   type(TYPE_NAME), intent(in)  :: this
   character(len=36)            :: str
   str = this%data%id
  end function id_

  function name_(this) result(str)
   type(TYPE_NAME), intent(in) :: this
   character(len=len_trim(this%data%name)) :: str
   str = trim(this%data%name)
  end function name_

  subroutine tag_new_object(this)
   type(TYPE_NAME), intent(inout)  :: this

   logical, external  :: print_debug_object_info
   external           :: get_uuid

    call get_uuid(this%data%id)
    if (print_debug_object_info()) then
      print *, '--> allocated ' // trim(this%data%name) // " " // id(this)
    endif

   end subroutine tag_new_object

!=============================================================
