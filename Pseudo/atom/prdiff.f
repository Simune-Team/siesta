C
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
 9000 format(/' &d total energy differences in series',
     $      //,' &d',2x,9i9)
 9010 format(' &d',1x,i2,1x,9f9.4)
      write(6,'(/,a,/)') '*----- End of series ----* spdfg &d&v'
c
      return
c
      end
