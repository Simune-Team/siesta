module external_interfaces
!
! The user must provide external subroutines
! with these interfaces
!
interface
   ! Called to terminate the program, optionally
   ! printing a message
      subroutine die(str)
      character(len=*), intent(in), optional :: str
      end subroutine die
end interface

end module external_interfaces
