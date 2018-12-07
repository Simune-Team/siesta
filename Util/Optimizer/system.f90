module system_m

  public :: system

  CONTAINS

!! Until we set F2008 as the reference, this is needed for F2003

subroutine system(str,stat)
  use iso_c_binding, only: C_CHAR, C_NULL_CHAR, C_INT

  character(len=*), intent(in)   :: str
  integer, intent(out), optional :: stat
   interface
      integer(c_int) function c_system(string) bind(C,name="system")
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char) :: string(*)
      end function c_system
   end interface

   integer(C_INT) :: return_value

   return_value = c_system(str // C_NULL_CHAR)
   if (present(stat)) then
      stat = return_value
   end if

 end subroutine system

end module system_m
