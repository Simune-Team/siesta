subroutine die(str)
character(len=*), intent(in) :: str
 write(0,*) trim(str)
 stop
end subroutine die

