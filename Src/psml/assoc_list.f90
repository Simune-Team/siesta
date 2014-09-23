module assoc_list
! First version, with fixed (initial) length,
! and fixed-length fields.
! Alberto Garcia, Sept 2014
!
!-----------------------------------------------------------
type, public :: assoc_list_t
  private
  integer                              :: nslots
  integer                              :: nitems
  character(len=50), allocatable       :: key(:)
  character(len=50), allocatable       :: value(:)
end type assoc_list_t

public :: assoc_list_init
public :: assoc_list_insert
public :: assoc_list_nitems
public :: assoc_list_get_key
public :: assoc_list_get_value

CONTAINS

subroutine assoc_list_init(a,n,stat)
type(assoc_list_t), intent(inout) :: a
integer, intent(in)               :: n
integer, intent(out)              :: stat

  if (allocated(a%key)) then
     deallocate(a%key)
  endif
  if (allocated(a%value)) then
     deallocate(a%value)
  endif
  a%nslots = n
  a%nitems = 0
  allocate(a%key(n),a%value(n),stat=stat)

end subroutine assoc_list_init

subroutine assoc_list_insert(a,key,value,stat)
type(assoc_list_t), intent(inout) :: a
character(len=*), intent(in)      :: key, value
integer, intent(out)              :: stat

integer :: i

! Replace if key exists already
do i = 1, a%nitems
   if (a%key(i) == key) then
      a%value(i) = value
      stat = 0
      return
   endif
enddo
!
! Add at the end
a%nitems = a%nitems + 1
if (a%nitems > a%nslots)  then
   stat = -1
   return
endif
i = a%nitems
a%key(i) = key
a%value(i) = value
stat = 0

end subroutine assoc_list_insert

function assoc_list_nitems(a) result(n)
type(assoc_list_t), intent(in)    :: a
integer                           :: n

n = a%nitems
end function assoc_list_nitems

subroutine assoc_list_get_key(a,i,key,stat)
type(assoc_list_t), intent(in)    :: a
integer, intent(in)               :: i
character(len=*), intent(out)     :: key
integer, intent(out)              :: stat

if (i > a%nitems) then
   stat = -1
   return
endif
key = a%key(i)
stat = 0
end subroutine assoc_list_get_key

subroutine assoc_list_get_value(a,key,value,stat)
type(assoc_list_t), intent(in) :: a
character(len=*), intent(in)      :: key
character(len=*), intent(out)     :: value
integer, intent(out)              :: stat

do i = 1, a%nitems
   if (a%key(i) == key) then
      value = a%value(i) 
      stat = 0
      return
   endif
enddo
stat = -1
end subroutine assoc_list_get_value

end module assoc_list
