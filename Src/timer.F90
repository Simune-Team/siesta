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

!!!!!!!!!!  AG: This should go in arch.make:  #define __PPC970__
!! Code by Rogeli Grima, BSC, specific for MareNostrum
!!
#ifdef __PPC970__

      MODULE tiempo
      integer(4)          :: maxcont = 20
      integer(4), pointer :: neventos(:)
      integer(8), pointer :: timer(:,:)
      real(8),    pointer :: abs_time(:)
      real(8),  parameter :: TRES=69.8356e-9
      END MODULE tiempo
 
      subroutine inicializar_tiempo( mxcnt )
      use tiempo
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(4), optional :: mxcnt
!--------------------------------------------------------------- Local Variables
      integer(4) :: ii, counter, rate, maxCounter
      integer(8) :: walltime
!------------------------------------------------------------------------- BEGIN
      if (present(mxcnt)) maxcont = mxcnt
      nullify( timer, abs_time, neventos )
      allocate( timer(2,maxcont), abs_time(maxcont), neventos(maxcont) )

      do ii= 1, maxcont
        timer(1,ii)  = 0
        timer(2,ii)  = 0
        neventos(ii) = 0.0
        abs_time(ii) = 0.0
      enddo
!--------------------------------------------------------------------------- END
      end subroutine inicializar_tiempo

      subroutine medirTiempo( job, cont )
      use tiempo
      use sys, only : die
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(4) :: job, cont
!--------------------------------------------------------------- Local Variables
      integer(8) :: itime, walltime

!------------------------------------------------------------------------- BEGIN
      if (cont < 1 .OR. cont > maxcont) then
        call die( 'Wrong cont value' )
      endif

      IF (job == 1) THEN
      	timer(1,cont) = walltime( )
      ELSE IF (job == 2) THEN
      	timer(2,cont)  = walltime( )
        itime          = timer(2,cont) - timer(1,cont)
        abs_time(cont) = abs_time(cont) + itime*TRES
        neventos(cont) = neventos(cont) + 1
      ELSE IF (job == 3) THEN
      	timer(2,cont)  = walltime( )
        itime          = timer(2,cont) - timer(1,cont)
        abs_time(cont) = abs_time(cont) + itime*TRES
        neventos(cont) = neventos(cont) + 1
#ifdef DEBUG
        write(23,*) '      iter:', neventos(cont),
     &    'Tiempo:', itime*TRES, ' total:', abs_time(cont)
#endif
      ELSE IF (job == 0) THEN
        timer(1,cont)  = walltime( )
        abs_time(cont) = 0.0
        neventos(cont) = 0.0
      ELSE
        call die( 'Wrong job value' )
      ENDIF
!--------------------------------------------------------------------------- END
      end subroutine medirTiempo

      subroutine showTiempo( uf, cont, mess )
      use tiempo
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(4)         :: uf, cont
      character(*)       :: mess
!--------------------------------------------------------------- Local Variables
      real(8)            :: tt
!------------------------------------------------------------------------- BEGIN
      tt = abs_time(cont)

      if (neventos(cont).eq.0) then
        write(uf,*) mess,  tt, neventos(cont), tt
      else
        write(uf,*) mess,  tt, neventos(cont), tt/neventos(cont)
      endif
!--------------------------------------------------------------------------- END
      end subroutine showTiempo
#else
      MODULE tiempo
      implicit none
      integer(4)          :: maxcont = 20
      integer(8), pointer :: timer(:,:)
      real(8),    pointer :: abs_time(:), neventos(:)
      integer(8)          :: current_time
      integer(8)          :: ref_time
      END MODULE tiempo

      subroutine inicializar_tiempo( mxcnt )
      use tiempo
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(4), optional :: mxcnt
!--------------------------------------------------------------- Local Variables
      integer(4):: ii, counter, rate, maxCounter
      real(8)   :: tt
!------------------------------------------------------------------------- BEGIN
      if (present(mxcnt)) maxcont = mxcnt
      nullify( timer, abs_time, neventos )
      allocate( timer(2,maxcont), abs_time(maxcont), neventos(maxcont) )

      do ii= 1, maxcont
        timer(1,ii)  = 0
        timer(2,ii)  = 0
        neventos(ii) = 0.0
        abs_time(ii) = 0.0
      enddo

      call system_clock( counter, rate, maxCounter )
      current_time = counter
      ref_time     = counter
      tt = REAL(maxCounter)/REAL(rate)
!--------------------------------------------------------------------------- END
      end subroutine inicializar_tiempo

      subroutine medirTiempo( job, cont )
      use tiempo
      use sys, only : die
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(4)           :: job, cont
!--------------------------------------------------------------- Local Variables
      integer(4)           :: rate, contMax, counter
      integer(8)           :: itime

!------------------------------------------------------------------------- BEGIN
      if (cont < 1 .OR. cont > maxcont) then
        call die( 'Wrong cont value' )
      endif

      call system_clock( counter, rate, contMax )
      if (counter .ge. ref_time) then
      	current_time = current_time + counter - ref_time
      else
        current_time = current_time + contMax - ref_time + counter + 1
      endif
      ref_time = counter

      IF (job == 1) THEN
      	timer(1,cont) = current_time
      ELSE IF (job == 2) THEN
      	timer(2,cont)  = current_time
        itime          = timer(2,cont) - timer(1,cont)
        abs_time(cont) = abs_time(cont) + (REAL(itime) / REAL(rate))
        neventos(cont) = neventos(cont) + 1.0
      ELSE IF (job == 3) THEN
      	timer(2,cont)  = current_time
        itime          = timer(2,cont) - timer(1,cont)
        abs_time(cont) = abs_time(cont) + (REAL(itime) / REAL(rate))
        neventos(cont) = neventos(cont) + 1.0
#ifdef DEBUG
        write(23,*) '      iter:', neventos(cont),  &
         'Tiempo:', (REAL(itime) / REAL(rate)),     &
         ' total:', abs_time(cont)
#endif
      ELSE IF (job == 0) THEN
        timer(1,cont)  = current_time
        abs_time(cont) = 0.0
        neventos(cont) = 0.0
      ELSE
        call die( 'Wrong job value' )
      ENDIF
!--------------------------------------------------------------------------- END
      end subroutine medirTiempo

      subroutine showTiempo( uf, cont, mess )
      use tiempo
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(4)         :: uf, cont
      character(*)       :: mess
!--------------------------------------------------------------- Local Variables
      real(8)            :: tt
      real(8), parameter :: ZERO=0.0
!------------------------------------------------------------------------- BEGIN
      tt = abs_time(cont)

      if (neventos(cont).eq.ZERO) then
        write(uf,*) mess,  tt, neventos(cont), tt
      else
        write(uf,*) mess,  tt, neventos(cont), tt/neventos(cont)
      endif
!--------------------------------------------------------------------------- END
      end subroutine showTiempo
#endif
