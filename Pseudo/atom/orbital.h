c------
c $Id: orbital.h,v 1.2 1997/05/22 17:32:23 wdpgaara Exp $
c
c $Log: orbital.h,v $
c Revision 1.2  1997/05/22 17:32:23  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
c
      integer norbmx
      parameter (norbmx=40)
c
c     norb:     Total number of orbitals.
c     ncp :     First valence orbital.
c
      integer norb, ncp
      character il(5)*1
c
      integer no(norbmx), lo(norbmx)
      double precision so(norbmx), zo(norbmx)
      logical down(norbmx)
c
      common /orbital/ so, zo
      common /orb_int/ norb, ncp, no, lo
      common /orb_char/ il
      common /orb_log/  down
      save /orbital/, /orb_int/, /orb_char/, /orb_log/
c------
