module m_wxml_dictionary

  use m_wxml_escape, only : check_Name

  implicit none

private
!
! A very rough implementation for now
! It uses fixed-length buffers for key/value pairs,
! and the maximum number of dictionary items is hardwired.

integer, parameter, private    :: MAX_ITEMS = 30

type, private :: dict_item
  character(len=1), pointer, dimension(:) :: key
  character(len=1), pointer, dimension(:) :: value
end type dict_item

type, public :: wxml_dictionary_t
private
      integer                               :: number_of_items ! = 0
      type(dict_item), dimension(MAX_ITEMS) :: items
end type wxml_dictionary_t

!
! Building procedures
!
public :: init_dict
public :: reset_dict
public :: add_item_to_dict
!
! Query and extraction procedures
!
public  :: len
interface len
   module procedure number_of_entries
end interface
public  :: number_of_entries
public  :: get_key, get_key_len
public  :: get_value, get_value_len
public  :: has_key
public  :: print_dict
!
interface get_value
   module procedure wxml_get_value
end interface

CONTAINS

!------------------------------------------------------
function number_of_entries(dict) result(n)
type(wxml_dictionary_t), intent(in)   :: dict
integer                          :: n

n = dict%number_of_items

end function number_of_entries

!------------------------------------------------------
function has_key(dict,key) result(found)
type(wxml_dictionary_t), intent(in)   :: dict
character(len=*), intent(in)     :: key
logical                          :: found

integer  :: n, i
found = .false.
n = dict%number_of_items
do  i = 1, dict%number_of_items
  if (size(dict%items(i)%key) == len(key)) then
    if (transfer(dict%items(i)%key,key) == key) then
      found = .true.
      exit
    endif
  endif
enddo
end function has_key

!------------------------------------------------------
subroutine wxml_get_value(dict,key,value,status)
type(wxml_dictionary_t), intent(in)            :: dict
character(len=*), intent(in)              :: key
character(len=*), intent(out)             :: value
integer, intent(out)                      :: status
!
integer  :: i

status = -1
do  i = 1, dict%number_of_items
  if (size(dict%items(i)%key) == len(key)) then
    if (transfer(dict%items(i)%key,key) == key) then
      value = transfer(dict%items(i)%value,value)
      status = 0
      exit
    endif
  endif
enddo

end subroutine wxml_get_value

function get_value_len(dict, key) result(value_len)
  type(wxml_dictionary_t), intent(in) :: dict
  character(len=*), intent(in) :: key
  integer :: value_len

  integer :: i

  value_len = 0
  do  i = 1, dict%number_of_items
    if (size(dict%items(i)%key) == len(key)) then
      if (transfer(dict%items(i)%key,key) == key) then
        value_len = size(dict%items(i)%value)
        exit
      endif
    endif
  enddo

end function get_value_len


!------------------------------------------------------
subroutine get_key(dict,i,key,status)
!
! Get the i'th key
!
type(wxml_dictionary_t), intent(in)            :: dict
integer, intent(in)                       :: i
character(len=*), intent(out)             :: key
integer, intent(out)                      :: status

if (i>0 .and. i<=dict%number_of_items) then
      key = transfer(dict%items(i)%key, key)
      status = 0
else
      key = ""
      status = -1
endif

end subroutine get_key

function get_key_len(dict, i) result(key_len)
  type(wxml_dictionary_t), intent(in) :: dict
  integer, intent(in) :: i
  integer :: key_len

if (i>0 .and. i<=dict%number_of_items) then
  key_len = size(dict%items(i)%key)
else
  key_len = 0
endif

end function get_key_len


subroutine add_item_to_dict(dict, key, value)

  type(wxml_dictionary_t), intent(inout) :: dict
  character(len=*), intent(in)          :: key
  character(len=*), intent(in)          :: value

  character(len=len(key)) :: check_key
  integer  :: n

  n = dict%number_of_items
  if (n == MAX_ITEMS) then
    write(unit=0,fmt=*) "Dictionary capacity exceeded !"
    RETURN
  endif

! keys may not have initial (or trailing; thus trim below) blanks:
!TOHW remove this check? shouldn't be passing blank-prefixed strings anyway.
  check_key=adjustl(key)
  if (.not.check_Name(trim(check_key))) then
    write(0,*) 'attribute name is invalid'
    call abort()
  endif

  n = n + 1
  allocate(dict%items(n)%key(len_trim(check_key)))
  dict%items(n)%key = transfer(trim(check_key),dict%items(n)%key)
  allocate(dict%items(n)%value(len(value)))
  dict%items(n)%value = transfer(value,dict%items(n)%value)

  dict%number_of_items = n

end subroutine add_item_to_dict

!------------------------------------------------------
subroutine init_dict(dict)
  type(wxml_dictionary_t), intent(out)   :: dict

  integer :: i
  
  do i = 1, MAX_ITEMS
    nullify(dict%items(i)%key)
    nullify(dict%items(i)%key)
  enddo

  dict % number_of_items = 0

end subroutine init_dict
  
!------------------------------------------------------
subroutine reset_dict(dict)
  type(wxml_dictionary_t), intent(inout)   :: dict
  
  integer :: i
  do i = 1, dict%number_of_items
    deallocate(dict%items(i)%key)
    deallocate(dict%items(i)%value)
  enddo

  dict%number_of_items = 0

end subroutine reset_dict

!------------------------------------------------------
subroutine print_dict(dict)
type(wxml_dictionary_t), intent(in)   :: dict

integer  :: i

do i = 1, dict%number_of_items
      print *, dict%items(i)%key, " = ", dict%items(i)%value
enddo

end subroutine print_dict

end module m_wxml_dictionary
