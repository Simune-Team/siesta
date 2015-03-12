! The PSML library calls a "die" routine when it encounters
! an error. This routine should take care of carrying out any
! needed cleaning and terminating the program.
! As the details would vary with the client program, each program
! has to provide its own.
! 
! This is an example implementation that could be used for serial programs.

subroutine die(str)
  character(len=*), intent(in), optional :: str

  if (present(str)) then
     write(0,"(a)") str
  endif
  STOP

end subroutine die
