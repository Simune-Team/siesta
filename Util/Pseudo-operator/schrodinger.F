
      module schrodinger_m

      use precision,        only: dp
      use atmparams,        only: nrmax

      implicit none

      CONTAINS

! ----------------------------------------------------------------------

      subroutine schro_eq( Zval, rofi, vps, ve, s, drdi,
     .                     nrc, l, a, b, nnodes, nprin,
     .                     e, g )

      implicit none

      integer 
     .  nrc,l,nnodes,nprin
      real(dp)
     .  Zval, rofi(*),vps(*),ve(*),s(nrc),drdi(*),a,b,e,g(*)
    
!     Internal variables

      real(dp)
     .  a2b4, h(nrmax), r2, vtot, rmax, dr,
     .  y(nrmax), dnrm, phi, dsq
      integer
     .  ir

      a2b4 = a * a * 0.25_dp

      do ir = 2, nrc
        g(ir) = 0.0_dp
        r2    = ( rofi(ir)**2 )
        vtot  = vps(ir) + ve(ir) + dble(l*(l+1))/r2
        h(ir) = vtot * s(ir) + a2b4
      enddo
      h(1) = h(2)
      g(1) = 0.0_dp

      e  = -( (zval/dble(nprin))**2 ) ! Guess for eigenvalue
      dr = -1.0d6                     ! Requested logarithmic derivative at Rmax
      rmax = rofi(nrc)
      call egofv( h, s, nrc, e, g, y, l, zval, a, b, rmax,
     .            nprin, nnodes, dr )

      do ir = 2, nrc
        phi   = g(ir)
        dsq   = sqrt( drdi(ir) )
        phi   = phi * dsq
        g(ir) = phi
      enddo
      g(1) = 0.0_dp
      dnrm = 0.0_dp

      do ir = 2, nrc
        phi  = g(ir)
        dnrm = dnrm + phi * phi * drdi(ir)
      enddo
      dnrm=sqrt(dnrm)

      do ir = 2, nrc
        g(ir) = g(ir) / dnrm
      enddo

      end subroutine schro_eq

! ----------------------------------------------------------------------
      subroutine energ_deriv( a, r, psi, vps,
     .                        ve, drdi, nrc, l, el, psidev, nrval )

      implicit none

!     This routine calculate the energy derivative of
!     a given wavefunction.
!     The routine solve and inhomogeneus version of
!     Schrodinger eqn.
!     It is not an optimized algorithm!!!!!!!!!!!!!!!!!
!     Written by Daniel Sanchez-Portal, July 1999

      integer  l, nrmin, nrval, ir, nrc

      real(dp)  r(nrval),psi(nrval),psidev(nrval),
     .  el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval),
     .  hi, dnrm, cons, a, ortog, dnrm2

      parameter(nrmin=1)
      nrc = min( nrc, nrval )

! Solving the inhomogeneus Schrodinger equation
      do ir = 2, nrc
        hi    = vps(ir) + ve(ir) + l*(l+1)/r(ir)**2-el
        hi    = hi * ( drdi(ir)**2 )
        hi    = hi + 0.25_dp * a**2
        h(ir) = hi
      enddo
      h(1) = h(2)

      cons = psi(nrmin+1) / ( vps(nrmin+1) + ve(nrmin+1) - el )
      cons = cons / r(nrmin+1)**(l+1)
      g(1) = 0.0_dp
      do ir = 1, nrmin+1
        g(ir) = cons * (r(ir)**(l+1))/sqrt(drdi(ir))
      enddo

      do ir = nrmin + 2, nrc
        hi = -( (psi(ir)+10.0_dp*psi(ir-1)
     .         + psi(ir-2))/12.0_dp )

        hi = hi + (10.0_dp*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2))/12.0_dp
        hi = hi + 2.0_dp * g(ir-1)-g(ir-2)

        g(ir) = hi / (1.0_dp-h(ir)/12.0_dp)
      enddo

! Orthogonalize the energy derivative to the original wavefunction
! and normalize
      dnrm2 = 0.0_dp
      ortog = 0.0_dp
      do ir = 1, nrc
        g(ir) = g(ir) * sqrt( drdi(ir) )
        dnrm2 = dnrm2 + drdi(ir) * (psi(ir)**2)
        ortog = ortog + drdi(ir) * g(ir) * psi(ir)
      enddo
      dnrm = 0.0_dp
      do ir = 1, nrc
        g(ir) = g(ir) - ortog * psi(ir) / dnrm2
        dnrm  = dnrm + drdi(ir) * (g(ir)**2)
      enddo
      dnrm = sqrt(dnrm)
      do ir = 1, nrc
        psidev(ir) = g(ir) / dnrm
      enddo

      end subroutine energ_deriv

! ----------------------------------------------------------------------
      subroutine rphi_vs_e( a, b, r, vps,
     .                      ve, nrval, l, el, rphi, rmax )

!      use basis_specs, only: restricted_grid

!   Calculate the atomic
!   radial wavefunction of the pseudopotential Vps, with angular
!   momentum  l, and energy el, inside r<Rmax
!   The Schrodinger equation is solved using a simple Numerov
!   scheme. Rmax should not be taken too big.
!   D. Sanchez-Portal, July 1999.


      real(dp) a, b
      integer nrval
      real(dp) r(nrval),
     .  el,vps(nrval),g(nrmax),drdi(nrmax),h(nrmax),ve(nrval),
     .  rphi(nrval), rmax, dnrm

      real(dp) big, dexpa, ab, hi
      parameter (big=1.0d6)
      integer  l, nrc, jr, ir


      dexpa = exp(a)
      ab    = a * b

      do ir = 1, nrval
        drdi(ir) = ab
        ab       = dexpa * ab
      enddo

      do ir = 2, nrval
        hi    = vps(ir) + ve(ir) + dble(l*(l+1))/r(ir)**2-el
        hi    = hi * ( drdi(ir)**2 )
        hi    = hi + 0.25_dp * a**2
        h(ir) = hi
      enddo
      h(1) = h(2)

      g(1) = 0.0_dp
      g(2) = 1.0_dp
      nrc = nint(log(rmax/b+1.0_dp)/a)+1
      nrc = min(nrc,nrval)
!      if (restricted_grid) nrc=nrc+1-mod(nrc,2)

      do ir = 3, nrc
        hi    = ( 10.0_dp*h(ir-1)*g(ir-1)+h(ir-2)*g(ir-2) )/12.0_dp
        hi    = hi + 2.0_dp * g(ir-1)-g(ir-2)
        g(ir) = hi / ( 1.0_dp-h(ir)/12.0_dp )

        if( abs(g(ir)) .gt. big ) then
           dnrm = 0.0_dp
           do jr = 1, ir
             dnrm = dnrm + drdi(jr) * ( g(jr) * sqrt(drdi(jr)) )**2
           enddo
           dnrm = sqrt(dnrm)
           do jr = 1, ir
             g(jr) = g(jr) / dnrm
           enddo
        endif
      enddo

! Normalize the wavefunction
      dnrm = 0.0_dp
      do ir = 1, nrc
        g(ir) = g(ir) * sqrt(drdi(ir))
        dnrm  = dnrm + drdi(ir) * ( g(ir)**2 )
      enddo

      dnrm = sqrt( dnrm )
      do ir = 1, nrc
        rphi(ir) = g(ir) / dnrm
      enddo

      end subroutine rphi_vs_e

! ----------------------------------------------------------------------

      end module schrodinger_m
