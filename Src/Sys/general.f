c $Id: general.f,v 1.8 1999/02/26 10:29:00 wdpgaara Exp $

      include 'cdiag_general.f'
      include 'eispack.f'
      include 'poison_general.f'
      include 'cft.f'
      include 'rdiag_general.f'

      SUBROUTINE CPUTIM (TIME)
      double precision time
c
c     Stub routine for systems without timing resources...
c
      TIME = 0.D0

      END

