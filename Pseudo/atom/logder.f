c
      subroutine logder(flag)
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
      include 'ode_blk.h'
      include 'ode_path.h'
c
c     Computes the logarithmic derivative (LD) as a function of energy.
c     The all-electron and pseudopotential results can be compared
c     to assess the transferability of the pseudopotential.
c
c     Note:  Non-relativistic wave functions are used.
c
c     Alberto Garcia, April 27, 1991
!     A. Garcia, May 2014
c
c-----
c     Number of energy values at which the LD is computed and
c     half energy range around the eigenvalue (~2 rydberg)
c
      integer npt_energ
      double precision half_range
      parameter (npt_energ=51,half_range=2.d0)
      double precision eps
      parameter (eps=1.d-7)
c
c     Flag: 'AE' for all-electron
c           'PS' for pseudo
c
      character flag*2
c
      double precision y(2)
c
      integer i, j, jr0, k, lp, llp, nok, nbad
      integer l, llo, lhi
      double precision emin, emax, eigv, step, r0, ld, h1, hmin
      character filename*9
c
      external ode, rkqc
c
c    Ode integration operational parameters
c
      h1 = 0.1d0
      hmin = 0.d0
      kmax = 200
      dxsav = 0.1d0
c
      llo = 10
      lhi = 0
      do i = ncp, norb
         if (lo(i) .gt. lhi) lhi = lo(i)
         if (lo(i) .lt. llo) llo = lo(i)
      enddo

!     loop over angular momenta in the valence shell

      do l = llo, lhi
!        Find the actual range of eigenvalues for this l

         write(filename,"(a2,a2,i1)") flag, 'EV', l
         open(unit=3,file=filename,form='formatted',status='unknown')

         emin = 1.0d10
         emax = -1.0d10
         do i = ncp, norb
            if (lo(i) .eq. l) then
               write(3,"(f12.6,f4.1)") ev(i), 0.0d0
               emin = min(emin,ev(i))
               emax = max(emax,ev(i))
            endif
         enddo
         close(3)

         emin = emin - half_range
         emax = emax + half_range
         step = (emax - emin) / (npt_energ - 1)
c
         lp = l + 1
         llp = l*lp
c
         write(filename,9900) flag, 'LOGD', l
 9900    format(a2,a4,i1)
         open(unit=3,file=filename,form='formatted',status='unknown')
c
c        Radius at which the LD will be calculated. Fix it at
c        a grid point. (Numerical Recipes, p.90)
c
         r0 = logder_radius
         call locate(r,nr,r0,jr0)
         r0 = r(jr0)

!        Use only down potential (Coulomb for AE, scalar rel for PS)
         do 50 j = 2, nr
            v(j) = viod(lp,j)/r(j) + vid(j) + llp/r(j)**2
 50      continue
c
        do 30 k = 0, npt_energ - 1
c
          energ = emin + step * k
c
c         Call the integration routine. Adaptive stepsize.
c
c         y(1): wavefunction   ( = r R )
c         y(2): first derivative
c
c         Initial values (at r(2), near r=0 and so R ~ r^l )
c
          y(1) = r(2) ** lp
          y(2) = lp * r(2)**(lp-1)
!
!         Integrate up to r0
!
          call odeint(y,2,r(2),r0,eps,h1,hmin,nok,nbad,ode,rkqc)
c
          ld = y(2) / y(1)

          write(3,8000) energ, ld

 8000     format(1x,f12.6,3x,f12.6)

 30       enddo
          close(3)

       enddo

      end

