!--------------------------------------------------------------
!
! Basic functionality for reference-counted data structures
!
!   This file has to be included after the global declarations
!   for extra functionality, just where a "contains" statement
!   would go
!
!   Parameter:  "TYPE_NAME" has to be set to the name of the type
!--------------------------------------------------------------

! Acknowledgements of inspiration from the community:
! 
!    Basic templating via include files as in the FLIBS project
!    by Arjen Markus
!
!    Implementation of the classic reference-counting paradigm
!    in the PyF95 project.

  PRIVATE

  public  :: TYPE_NAME
  public  :: init, delete, assignment(=), refcount, id

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

  ! Stand-alone routine which must be provided
  interface 
      subroutine die(str)
      character(len=*), intent(in), optional  :: str
      end subroutine die
  end interface

CONTAINS

! -- Main structural features
!
   subroutine init_(this)

   ! Initializes new storage

     type (TYPE_NAME), intent(inout) :: this

     integer :: error

     ! First, remove the current reference
     call delete(this)

     ! Allocate fresh storage
     allocate(this%data, stat=error)
     if (error /= 0) then
       call die("Error allocating data structure")
     endif

     ! Set the initial reference count
     this%data%refCount = 1

  end subroutine init_

  subroutine delete_(this)

    ! Removes the current reference, possibly
    ! deallocating storage

    type (TYPE_NAME), intent(inout) :: this

    integer :: error
!    logical, external  :: print_debug_object_info

    if (.not. associated(this%data)) return

    this%data%refCount = this%data%refCount - 1

    if (this%data%refCount == 0) then

      ! Safe to delete the data now
      ! Use the routine provided for this specific
      ! type to clean up any internal structure

      call delete_Data(this%data)

!      if (print_debug_object_info()) then
!        print *, "--> deallocated " // id(this) // " " // trim(this%data%name) 
!      endif

      ! Deallocate the currently referenced storage

      deallocate(this%data, stat=error)
      if (error /= 0) then
         call die("Error in deallocation")
      endif
    endif

    this%data => null()

  end subroutine delete_


  subroutine assign_(this, other)

  ! Make "this" reference the same data  as "other".
  ! No copying of data is involved, simply an increment of the
  ! reference counter.

  ! IMPORTANT NOTE: Assignment must take the form of a subroutine, and
  ! not of a function, since the "inout" intent is essential. One has
  ! to clean up "this" before making it point to the same place as
  ! "other".  In a function, the intrinsic "out" intent for "this"
  ! will destroy any prior information.

    type (TYPE_NAME), intent(inout) :: this
    type (TYPE_NAME), intent(in) :: other

    if (.not. associated(other%data)) then
     call die("Assignment of non-initialized object in " // trim(mod_name))
    endif

    ! Delete to remove the current reference of "this"

    call delete(this)

    ! Establish the new reference and increment the reference counter.

    this%data => other%data
    this%data%refcount = this%data%refcount+1

  end subroutine assign_

  function refcount_(this) result(count)
   type(TYPE_NAME), intent(in)  :: this
   integer  :: count
   count = this%data%refCount
  end function refcount_

  ! The remaining procedures are not essential

  function id_(this) result(str)
   type(TYPE_NAME), intent(in)  :: this
   character(len=36)            :: str
   str = this%data%id
  end function id_

  subroutine tag_new_object(this)
   type(TYPE_NAME), intent(inout)  :: this

!   logical, external  :: print_debug_object_info
!   external           :: get_uuid

!    call get_uuid(this%data%id)
!    if (print_debug_object_info()) then
!      print *, '--> allocated ' // trim(this%data%name) // " " // id(this)
!    endif

   end subroutine tag_new_object

!=============================================================
