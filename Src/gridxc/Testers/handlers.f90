! Handlers
! ------------------------------------
!
! Example of routine die
! In your program, depending on your own usage details,
! you might need something more sophisticated
!
subroutine die(str)
 character(len=*), intent(in) :: str

 write(*,*) trim(str)
 stop
end subroutine die
!
! Handlers for alloc
! In your program, depending on your own usage details,
! you might need something more sophisticated
!
subroutine alloc_memory_event(bytes,name)
integer, intent(in) :: bytes
character(len=*), intent(in) :: name
!! write(*,*) "alloc: allocated ", bytes, "bytes for "//trim(name)
end subroutine alloc_memory_event

subroutine alloc_error_report(name,code)
character(len=*), intent(in) :: name
integer, intent(in) :: code
!! write(*,*) "alloc error: "//trim(name)
end subroutine alloc_error_report
!
! Timer interface: called at specific events
! 
         subroutine gridxc_timer_start(str)
         character(len=*), intent(in)  :: str
         end subroutine gridxc_timer_start

         subroutine gridxc_timer_stop(str)
         character(len=*), intent(in)  :: str
         end subroutine gridxc_timer_stop
!
! Include this here to dispatch internal timing requests
! (as they might be needed by the code for load-balancing)

!=============================================================================
! This module should be merged into whatever gridxc options
! module is finally used, and made accessible to the user
! via some 'gridxc_init' call
!
module timer_options
 logical, public :: use_walltime = .false.
 logical, public :: time_mpi_calls = .false.
end module timer_options

subroutine timer( prog, iOpt )

!  FINDS AND PRINTS THE CPU TIME SPENT IN DIFFERENT ROUTINES AND/OR
!   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=1 AT THE BEGINNING
!   OF EACH ROUTINE AND WITH IOPT=2 AT THE END OF IT.
!  ARGUMENTS:
!    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
!    IOPT: INPUT OPTION PARAMETER:
!      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
!      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
!      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
!      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
!  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
!  WRITTEN BY J.SOLER. JUL.2009

  use m_timer, only: timer_init    ! Initialize all times
  use m_timer, only: timer_start   ! Start counting time
  use m_timer, only: timer_stop    ! Stop counting time

! Arguments
  implicit none
  character(len=*),intent(in):: prog   ! Name of program to time
  integer,         intent(in):: iOpt   ! Action option

! Select action
  if (iOpt==0) then
     call timer_init()
  else if (iOpt==1) then
    call timer_start( prog )
  else if (iOpt==2) then
    call timer_stop( prog )
  else if (iOpt==3) then
   ! do nothing here
  else
    call die('timer: ERROR: invalid iOpt value')
  end if

end subroutine timer

subroutine timer_mpi( prog, iOpt )
 use timer_options, only : time_mpi_calls

  implicit none
  character(len=*),intent(in):: prog   ! Name of program to time
  integer,         intent(in):: iOpt   ! Action option

  if (time_mpi_calls) then
     call timer(prog,iOpt)
  endif

end subroutine timer_mpi

function use_walltime_in_timer()
  use timer_options, only: use_walltime

  logical :: use_walltime_in_timer
  use_walltime_in_timer = use_walltime
end function use_walltime_in_timer

