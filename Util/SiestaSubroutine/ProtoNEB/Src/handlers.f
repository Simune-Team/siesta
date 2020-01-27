
      subroutine alloc_error_report(str,code)

      character(len=*), intent(in)  :: str
      integer, intent(in)  :: code

      external :: die
      
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
!      use memory_log, only: memory_event

      integer, intent(in)           :: bytes
      character(len=*), intent(in)  :: name

      ! call memory_event(bytes,name)

      end subroutine alloc_memory_event

! The PSML library calls a "die" routine when it encounters an
! error. This routine should take care of carrying out any needed
! cleaning and terminating the program.  As the details would vary with
! the client program, each program has to provide its own.
! 

      subroutine psml_die(str)
      character(len=*), intent(in) :: str

      external :: die
      
      call die(str)

      end subroutine psml_die
