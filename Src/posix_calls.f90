module posix_calls

!!! See also
!!! https://stackoverflow.com/questions/30279228/is-there-an-alternative-to-getcwd-in-fortran-2003-2008
!!! for "Windows" code
  
  public :: system
  public :: chdir
  public :: getcwd

  private

contains

!! Until we set F2008 as the reference, this is needed for F2003

  !< Run a system command via `str` and possibly return the status value for the call
  subroutine system(str, stat)
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
    if ( present(stat) ) then
      stat = return_value
    end if

  end subroutine system

  !< Change directory in the execution
  subroutine chdir(path, stat)
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR, C_INT
    
    character(*) :: path
    integer, optional, intent(out) :: stat
    
    integer(C_INT) :: return_value
    
    interface
      ! int chdir(const char *path);
      integer(c_int) function c_chdir(path) bind(C,name="chdir")
        use iso_c_binding, only: c_char, c_int
        character(kind=c_char) :: path(*)
      end function c_chdir
    end interface
    
    return_value = c_chdir(path//c_null_char)
    
    if ( present(stat) ) then
      stat = return_value
    end if
    
  end subroutine chdir

!!$  !< Create a new directory, default to mode 0755
!!$  subroutine mkdir(path, mode, stat)
!!$    use iso_c_binding, only: C_CHAR, C_NULL_CHAR, C_INT
!!$    
!!$    character(*) :: path
!!$    integer, optional, intent(in) :: mode
!!$    integer, optional, intent(out) :: stat
!!$
!!$    integer(C_INT) :: use_mode, return_value
!!$
!!$    interface
!!$      ! int mkdir(const char *pathname, mode_t mode);
!!$      integer(c_int) function c_mkdir(path, mode) bind(C,name="mkdir")
!!$        use iso_c_binding, only: c_char, c_int
!!$        character(kind=c_char) :: path(*)
!!$        character(kind=c_int), value :: mode
!!$      end function c_mkdir
!!$    end interface
!!$
!!$    ! Determine directory mode
!!$    if ( present(mode) ) then
!!$      use_mode = int(o'775', c_int)
!!$    else
!!$      use_mode = mode
!!$    end if
!!$    return_value = c_mkdir(path//c_null_char, use_mode)
!!$    
!!$    if ( present(stat) ) then
!!$      stat = return_value
!!$    end if
!!$
!!$  end subroutine mkdir

  !< Query the current directory
  !<
  !< If `path` is too short to hold the result stat will be returned with value `1`.
  !< Otherwise `stat` will be `0`.
  !< You should prefer to use `stat`.
  subroutine getcwd(path, stat)
    use iso_c_binding, only: C_CHAR, C_INT, C_ASSOCIATED, C_PTR, C_SIZE_T
    
    character(len=*), intent(out)  :: path
    integer, optional, intent(out) :: stat
    
    interface
      ! char * getcwd(char *buf, size_t size);
      function c_getcwd(str,str_size) result(res) bind(C, name="getcwd")
        use iso_c_binding, only: C_CHAR, C_PTR, C_SIZE_T
        character(kind=C_CHAR), intent(out) :: str(*)
        integer(kind=C_SIZE_T), intent(in), VALUE  :: str_size
        type(C_PTR) :: res
      end function c_getcwd
    end interface
    
    integer(C_INT) :: return_value
    
    integer :: null_loc
    integer(kind=C_SIZE_T) :: str_size
    type(C_PTR) :: p_result
    
    str_size = len(path)
    p_result = c_getcwd(path, str_size)
    
    if ( c_associated(p_result) ) then
      ! Remove stuff past the null-character
      null_loc = index( path, char(0) )
      path = path(1:null_loc-1)
      return_value = 0
    else
      path = ""
      return_value = 1
    endif

    if ( present(stat) ) then
      stat = return_value
    endif

  end subroutine getcwd
  
end module posix_calls
