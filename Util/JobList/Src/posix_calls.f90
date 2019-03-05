module posix_calls

  public :: system
  public :: chdir
  public :: getcwd
  
  private

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
 
!-----------------------------------------------------
  subroutine chdir(path, stat)
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR, C_INT


    character(*) :: path
    integer, optional, intent(out) :: stat

    integer(C_INT) :: return_value

  interface
    integer(c_int) function c_chdir(path) bind(C,name="chdir")
      use iso_c_binding, only: c_char, c_int
      character(kind=c_char) :: path(*)
    end function c_chdir
  end interface

  return_value =  c_chdir(path//c_null_char)

  if (present(stat)) then
     stat = return_value
  endif
  end subroutine chdir
!--------------------------------------------------

  ! If 'path' is too short to hold the result, stat will be "1".
  ! Otherwise, stat=0
  ! You *should* use the stat argument!
  
  subroutine getcwd(path, stat)
    use iso_c_binding, only: C_CHAR, C_INT, C_ASSOCIATED, C_PTR, C_SIZE_T

    character(len=*), intent(out)  :: path
    integer, optional, intent(out) :: stat

interface
   !      char * getcwd(char *buf, size_t size);

   function c_getcwd(str,str_size) result(res) bind(C, name="getcwd")
     use ISO_C_Binding, only: C_CHAR, C_PTR, C_SIZE_T
     character(kind=C_CHAR),intent(out) :: str(*)
     integer(kind=C_SIZE_T),intent(in), VALUE  :: str_size
     type(C_PTR) :: res
    end function
  end interface

  integer(C_INT) :: return_value
  
  integer :: null_loc
  integer(kind=C_SIZE_T) :: str_size
  type(c_PTR) :: p_result
  
  str_size = len(path)
  p_result = c_getcwd(path,str_size)

  if (c_associated(p_result)) then
     ! Remove stuff past the null-character
     null_loc = index( path, char(0) )
     path = path(1:null_loc-1)
     return_value = 0
  else
     path = ""
     return_value = 1
  endif
  
  if (present(stat)) then
     stat = return_value
  endif
  
end subroutine getcwd

!!! See also
!!! https://stackoverflow.com/questions/30279228/is-there-an-alternative-to-getcwd-in-fortran-2003-2008
!!! for "Windows" code

  
end module posix_calls
