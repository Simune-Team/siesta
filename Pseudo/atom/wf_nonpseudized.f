c
      subroutine wf_nonpseudized(i,ar,br)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
c     Solves the wave equation for orbital i
c
      double precision zero, ai, pnine
      parameter (zero=0.d0,ai=2*137.0360411D0,pnine=0.9d0)
c
      integer i, nextr
      double precision rextr, rzero
c
      double precision ar(*), br(*)
c
      integer iflag, ist, j, lp, llp, ka
      double precision arp, arpm
      double precision v(nrmax)
c
      lp = lo(i) + 1
      llp = lo(i)*lp
c
      do 10 j = 1, nr
         ar(j) = 0.d0
   10 continue
      if (down(i)) then
         do 20 j = 2, nr
            v(j) = viod(lp,j)/r(j) + vid(j)
   20    continue
      else
         do 30 j = 2, nr
            v(j) = viou(lp,j)/r(j) + viu(j)
   30    continue
      end if
c
      if ( .not. relativistic) then
c
c           Add 'centrifugal term'
c
         do 40 j = 2, nr
            v(j) = v(j) + llp/r(j)**2
   40    continue
      end if
c
      if ( .not. relativistic) then
         call difnrl(0,i,v,ar,br,no(i),lo(i),so(i),ev(i),iflag)
      else
         call difrel(0,i,v,ar,br,no(i),lo(i),so(i),ev(i))
      end if
c
c     Plot and make the wavefunction 'upright'
c
!      ist = nint(sign(1.d0,ar(nr-85)))
!      call potrw(ar,r,nr-85,lo(i),1,ist,rc(lp))
c
!      do 50 j = 1, nr
!         ar(j) = ar(j)*ist
!         br(j) = br(j)*ist
!   50 continue

      return
c
      end


