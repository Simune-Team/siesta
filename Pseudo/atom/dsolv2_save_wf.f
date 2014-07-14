c
      subroutine dsolv2_save_wf(iter,id,i,nn,wfn)
c
c     This version serves to get the "major" "test" pseudo-wave-functions, 
c     that is, the solutions for *all* the valence shells with the
c     "down" pseudopotentials.
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
c
c  dsolv2 finds the (non) relativistic wave function using
c  difnrl to integrate the Schroedinger equation or
c  difrel to integrate the Dirac equation.
c  The energy level from the previous iteration are used
c  as initial guesses, and it must therefore be reasonably
c  accurate.
c
C     .. Parameters ..
      double precision zero, smev
      parameter (zero=0.D0,smev=1.D-4)
C     ..
C     .. Scalar Arguments ..
      integer iter, i
      character id*1
C     ..
c
c     Watch out for nn: it is not always equal to 'no'. 
c     In particular, some nn(i) could be zero.
c     The same is true for id...
c
      integer nn(*)
      double precision  wfn(nrmax)

C     ..
C     .. Local Scalars ..
      integer iflag, j, llp, lp
C     ..
C     .. External Subroutines ..
      external difnrl, difrel, orban
      logical leqi
      external leqi
C     ..
      double precision ar(nrmax), br(nrmax), v(nrmax)
C     ..

      if (ev(i) .ge. 0.0D0) ev(i) = -smev
c
c  Set up the potential, set the wave function array to zero-ar.
c
         lp = lo(i) + 1
         llp = lo(i)*lp
         do 40 j = 1, nr
            ar(j) = zero
   40    continue
         if (down(i)) then
            do 50 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j)
   50       continue
         else
            do 60 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j)
   60       continue
         end if
         if (.not. leqi(id,'r')) then
            do 70 j = 2, nr
               v(j) = v(j) + llp/r(j)**2
   70       continue
         end if
c
c  Call the integration routine.
c
         if (.not. leqi(id,'r')) then
            call difnrl(iter,i,v,ar,br,nn(i),lo(i),so(i),ev(i),iflag)
         else
            call difrel(iter,i,v,ar,br,nn(i),lo(i),so(i),ev(i))
         end if
c
         do 320 j = 1, nr
            wfn(j) = ar(j)
 320     continue
c
      return
c
      end
