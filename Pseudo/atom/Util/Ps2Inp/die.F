! To be used for serial programs, with a proper interface block

      subroutine die(str)
      character(len=*), intent(in), optional :: str
      if (present(str)) then
         write(*,*) str
      endif
      STOP
      end subroutine die
