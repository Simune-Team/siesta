c
      subroutine potrw(vd,r,nr,k,kj,ist,rc)
c
c  Step size of 0.05 is adjustable as seen fit to give
c  a reasonalble plot.
c
C     .. Parameters ..
      double precision zero, pzf
      parameter (zero=0.D0,pzf=0.05D0)
C     ..
C     .. Scalar Arguments ..
      integer ist, k, kj, nr
      double precision rc
C     ..
C     .. Array Arguments ..
      double precision r(nr), vd(nr)
C     ..
C     .. Local Scalars ..
      double precision step
      integer j
      character filename*7
C     ..
      if (kj .eq. 0) then
         write(filename,9900) k
 9900 format('PSWFNR',i1)
      else
         write(filename,9910) k
 9910 format('AEWFNR',i1)
      endif
c
      open(unit=3,file=filename,form='formatted',status='unknown')
c
c     Write out r, the wavefunction, and rc (kludge to pass it to
c     the plotting program)
c
      step = zero
      do 10 j = 2, nr
         if (r(j) .ge. step) then
            write(3,9000) r(j), vd(j)*ist, rc
            step = step + pzf
         end if
   10 continue
 9000 format(1x,f7.4,3x,2f10.6)
c
      close(unit=3)
cag
c      if (kj .eq. 0) then
c         if (k .eq. 0) then
c            marker = 'wsp'
c         else if (k .eq. 1) then
c            marker = 'wpp'
c         else if (k .eq. 2) then
c            marker = 'wdp'
c         else
c            marker = 'wfp'
c         end if
c      else
c         if (k .eq. 0) then
c            marker = 'wst'
c         else if (k .eq. 1) then
c            marker = 'wpt'
c         else if (k .eq. 2) then
c            marker = 'wdt'
c         else
c            marker = 'wft'
c         end if
c      end if
c      write(3,9010) marker
c 9010 format(1x,'marker ',a3)
c
      return
c
      end
