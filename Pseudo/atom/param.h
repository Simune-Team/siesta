c------
c $Id: param.h,v 1.3 1997/05/22 17:32:23 wdpgaara Exp $
c
c $Log: param.h,v $
c Revision 1.3  1997/05/22 17:32:23  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.2  1992/04/14  18:03:55  alberto
c Added logder_radius.
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
c     General operating parameters
c
      character ispp*1, icorr*2, nameat*2
      character irel*3, nicore*4
      double precision rsh, zel, znuc, zsh
c
      integer indd(5), indu(5)
      double precision rc(5), rc_input(5), cfac, rcfac,
     &                 logder_radius
      integer scheme, ncore, job, ifcore
      logical normal, polarized, relativistic
c
      common /param/ scheme, job, ifcore, ncore,
     &               rsh, zel, znuc, zsh, rc, rc_input, cfac, rcfac,
     &               logder_radius, indu, indd
      common /par_char/ nameat, ispp, icorr, irel, nicore
      common /par_log/  normal, polarized, relativistic
      save /param/, /par_char/, /par_log/
c------
