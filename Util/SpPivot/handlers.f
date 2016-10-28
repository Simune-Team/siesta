!--------------------------------------------------
! Stand-alone 'die' routine for use by libraries and
! low-level modules.
!
! Each program using the module or library needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta:
!
!   - The use of a Siesta-specific 'mpi_siesta' module.
!   - The need to have the pxf functionality.
!   - The use of 'unit 6' as output.
!
! Routines using this functionality should include
! the following
!
!     interface
!      subroutine die(str)
!      character(len=*), intent(in)  :: str
!      end subroutine die
!     end interface
!
!------------------------------------------------------

      subroutine die(str)

      character(len=*), intent(in)  :: str

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)
      write(6,'(a,i4)') 'Stopping Program from Node: ', Node
      write(0,'(a,i4)') 'Stopping Program from Node: ', Node
         call pxfflush(6)
         call pxfflush(0)
      call pxfabort()

      end subroutine die

      subroutine timer(str,i)

      character(len=*), intent(in)  :: str
      integer,  intent(in)  :: i
      end subroutine timer

!--------------------------------------------------
! Stand-alone routine to capture error messages from
! the alloc module
!
! This functionality could be made more general,
! and use a uniform interface for all the utility
! modules developed in-house. (Let's say, call it
! 'error_report' with severity arguments, etc).
!
! Each program using the alloc module needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta:
!
!   - The tagging of the message by Node.
!   - The use of 'unit 6' as output and '0' as error.
!
! Routines using this functionality should include
! the following
!
!   subroutine alloc_error_report(str,code)
!     character(len=*), intent(in) :: str
!     integer, intent(in)          :: code
!   end subroutine alloc_error_report
!
!------------------------------------------------------

      subroutine alloc_error_report(str,code)

      character(len=*), intent(in)  :: str
      integer, intent(in)  :: code

      if (code == 0) then
        call die(str)
      else
         write(0,*) trim(str)
         write(6,*) trim(str)
      endif

      end subroutine alloc_error_report
!--------------------------------------------------
! Stand-alone routine to capture memory events from
! the alloc module
!
! Each program using the alloc module needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta
!
!   - The use of the memory_log module for reports
!
! Routines using this functionality should include
! the following
!
!   subroutine alloc_memory_event(str,code)
!     character(len=*), intent(in) :: str
!     integer, intent(in)          :: code
!   end subroutine alloc_memory_event
!
!------------------------------------------------------

      subroutine alloc_memory_event(bytes,name)

      integer, intent(in)           :: bytes
      character(len=*), intent(in)  :: name

      ! Do nothing

      end subroutine alloc_memory_event

! The PSML library calls a "die" routine when it encounters an
! error. This routine should take care of carrying out any needed
! cleaning and terminating the program.  As the details would vary with
! the client program, each program has to provide its own.
! 

      subroutine psml_die(str)
      character(len=*), intent(in) :: str

      write(0,"(a)") str
      STOP

      end subroutine psml_die
