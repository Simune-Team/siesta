c
      subroutine denplot
c
c  Prints the charge density
c     
      include 'radial.h'
      include 'charge.h'
c
      double precision pi
      parameter (pi=3.141592653589d0)
c
      double precision fx
      integer j
c
      open(unit=3,file='CHARGE',form='formatted',status='unknown')
      open(unit=4,file='RHO',form='formatted',status='unknown')
c
c     Write out r, cdu, cdd and cdc
c

      do 10 j = 2, nr
            write(3,9000) r(j), cdu(j), cdd(j), cdc(j)
            fx = 1.d0 / ( 4.d0 * pi * r(j)**2)
            write(4,9000) r(j), fx*cdu(j), fx*cdd(j), fx*cdc(j)
   10 continue
 9000 format(1x,f15.10,3x,3f18.5)
c
      close(unit=3)
c
      return
c
      end
