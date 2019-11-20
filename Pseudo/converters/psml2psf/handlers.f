!--------------------------------------------------
! Stand-alone 'die' routine for use by libraries and
! low-level modules.
!
! Each program using the module or library needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
!
! Routines using this functionality might include
! the following for extra checking by the compiler
!
!     interface
!      subroutine die(str)
!      character(len=*), intent(in)  :: str
!      end subroutine die
!     end interface
!
! but note that a simple "external" statement, or even
! nothing at all, would work. This is because we have
! removed the "optional" character of the 'str' argument.
!------------------------------------------------------

! This is a bare-bones version with the basic functionality

      subroutine die(str)

      character(len=*), intent(in) :: str

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)
      STOP

      end subroutine die
      subroutine psml_die(str)

      character(len=*), intent(in) :: str

      call die(str)

      end subroutine psml_die
