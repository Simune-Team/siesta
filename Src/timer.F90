! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
!===============================================================================
! This is now a wrapper for the functionality in module m_timer
!
! subroutine timer( prog, iOpt )
!
!  FINDS AND PRINTS THE CPU TIME SPENT IN DIFFERENT ROUTINES AND/OR
!   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=0 AT THE BEGINNING
!   OF EACH ROUTINE AND WITH IOPT=1 AT THE END OF IT.
!  ARGUMENTS:
!    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
!    IOPT: INPUT OPTION PARAMETER:
!      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
!      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
!      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
!      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
!  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
!  WRITTEN BY J.SOLER. JUL.2009
!===============================================================================

subroutine timer( prog, iOpt )

! Module procedures used
  use sys,     only: die           ! Termination routine
  use m_timer, only: timer_init    ! Initialize all times
  use m_timer, only: timer_start   ! Start counting time
  use m_timer, only: timer_stop    ! Stop counting time
  use m_timer, only: timer_report  ! Write all times

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
    call timer_report( prog, printNow=.true. )
  else
    call die('timer: ERROR: invalid iOpt value')
  end if

end subroutine timer

