module system_m

  public :: system

  CONTAINS

!! Until we set F2008 as the reference, this is needed for F2003

subroutine system(str)
  use iso_c_binding, only: C_CHAR, C_NULL_CHAR

 character(len=*), intent(in) :: str
 interface
    integer(c_int) function my_system(string) bind(C,name="system")
      use iso_c_binding, only: c_char, c_int
      character(kind=c_char) :: string(*)
    end function my_system
 end interface

   ! Currently we ignore the return value.
   integer :: ignored_return_value

   ignored_return_value = my_system(str // C_NULL_CHAR)

 end subroutine system
end module system_m
