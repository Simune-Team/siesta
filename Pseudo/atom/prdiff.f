C
c $Id: prdiff.f,v 1.2 1997/05/22 17:32:27 wdpgaara Exp $
c
c $Log: prdiff.f,v $
c Revision 1.2  1997/05/22 17:32:27  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:55  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine prdiff(nconf,econf)
c
c   Prints out the energy differences between
c   different atomic configurations.
c
c   njtj  ***  modifications  ***
c     econf is able to handle larger numbers
c     of configurations.
c   njtj  ***  modifications  ***
c
c
C     .. Scalar Arguments ..
      integer nconf
C     ..
C     .. Array Arguments ..
      double precision econf(100)
C     ..
C     .. Local Scalars ..
      integer i, j
C     ..
      write(6,9000) (i,i=1,nconf)
      do 10 i = 1, nconf
         write(6,9010) i, (econf(i)-econf(j),j=1,i)
   10 continue
 9000 format(/' total energy difference',//2x,9i9)
 9010 format(1x,i2,1x,9f9.4)
c
      return
c
      end
