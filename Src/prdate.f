C $Id: prdate.f,v 1.5 1999/02/26 10:06:11 wdpgaara Exp $

      SUBROUTINE PRDATE (PROG)

C PRINTS A HEADING WITH THE SYSTEM AND PROGRAM NAMES, DATE AND TIME.
C THIS VERSION WORKS IN MOST UNIX WORKSTATIONS.
C WRITTEN BY P.SERENA AND J.M.SOLER

      CHARACTER PROG*(*)
*     CHARACTER HOST*20, ALLDAT*60

*     CALL SYSTEM('hostname > /tmp/prdate')
*     CALL SYSTEM('date >> /tmp/prdate')
*     OPEN(33, FILE='/tmp/prdate', STATUS='OLD')
*     READ(33,'(A20)') HOST
*     READ(33,'(A60)') ALLDAT
*     CLOSE(33)
*     CALL SYSTEM('rm /tmp/prdate')
*     WRITE(6,'(/,4A,/,2A,/)')
*    .   'prdate: PROGRAM = ', PROG, ', SYSTEM = ', HOST,
*    .   'prdate: DATE = ', ALLDAT

      END
