C $Id: cray.f,v 1.8 1999/02/26 18:04:47 wdpgaara Exp $

      include 'cdiag_general.f'
      include 'cft.f'
      include 'eispack.f'
      include 'poison_general.f'
      include 'rdiag_general.f'

      SUBROUTINE CPUTIM (TIME)

C  RETURNS CPU TIME IN SECONDS SINCE PROGRAM START
C  WRITEN BY J.SOLER (JSOLER AT EMDUAM11)

      DOUBLE PRECISION TIME
*     REAL TIMES(2)

C  NEXT LINES FOR IBM SYSTEMS-370 (ASSEMBLE ROUTINE TASKTM REQUIRED)
*     CALL TASKTM (ITIME)
*     TIME = 1.D-4 * ITIME

C  NEXT LINES FOR IBM-3090
*     CALL CPUTIME (TIME,RCODE)
*     TIME = 1.D-6 * TIME

C  NEXT LINE FOR CRAY
      TIME = SECOND()

C  NEXT LINE FOR SUN OR DEC WORKSTATIONS
*     TIME = ETIME(TIMES)

C  NEXT LINE FOR IBM RS/6000 WORKSTATIONS
*     TIME = MCLOCK()*0.01D0

C  NEXT LINE FOR ANYTHING ELSE
cc      TIME = 0.D0

      END
