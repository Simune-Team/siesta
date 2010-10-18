c $Header:$
c**********************************************************************
c Armiento Mattsson am05 functional for exchange and correlation,
c  
c Spin polarization type: exchange index scales with the separate 
c spin densities; correlation index and prefactor both scales with
c the separate spin densities. (xscss)
c 
c Version: 7 (xscss, web)
c
c The latest version of this subroutine file is available at
c   http://dft.sandia.gov/functionals/AM05.html
c
c
c Usage:
c
c Below follows first a number of usage examples on how to calculate energy
c and potential in codes that use different schemes for calculating
c the potential, and therefore uses different input quantities. After the
c usage examples follows the main am05 routine.
c
c Examples for the following schemes are included:
c
c   am05trads : The traditional scheme used e.g. by PRB 33, 8800 (1986).
c   am05wbs : The White and Bird scheme [PRB 50, 4954 (1994)].
c   am05pgjs : The Pople, Gill, and Johnson scheme [CPL 199, 557 (1992)].
c   am05wbnums : The White and Bird scheme with numerical derivatives.
c
c
c Citation request: 
c
c When using this functional, please cite:
c "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005);
c  R. Armiento and A. E. Mattsson (unpublished)."
c (The first paper for the AM05 functional, the second for the
c spin-polarized version)
c
c
c License and copyright:
c
c (c) Rickard Armiento 2005-2008
c
c Permission is hereby granted, free of charge, to any person obtaining 
c a copy of this software and associated documentation files (the 
c "Software"), to deal in the Software without restriction, including 
c without limitation the rights to use, copy, modify, merge, publish, 
c distribute, sublicense, and/or sell copies of the Software, and to 
c permit persons to whom the Software is furnished to do so, subject 
c to the following conditions:
c
c The above copyright notice and this permission notice shall be 
c included in all copies or substantial portions of the Software.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
c OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
c NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
c HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
c WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
c FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
c OTHER DEALINGS IN THE SOFTWARE.
c
c**********************************************************************

c**********************************************************************
c saferecp
c Helper function for making divisions with built in cutoff 
c for very small denominators
c
c output
c   set   variable to assign the value of the reciprocal
c
c input
c   nom   nominator
c   denom denominator
c**********************************************************************

      subroutine saferecp(set, nom, denom)

c     ** Input parameters
      real*8 nom, denom

c     ** Output parameters
      real*8 set

      if(denom .ge. 1e-30) then
         set = nom/denom
      else if(nom .le. 1e-30) then
         set = 0.0d0
      else if(nom .ge. 1e-30) then
         set = 1d30
      endif

      return
      end

c**********************************************************************
c am05wbs
c Usage example for White and Bird scheme [PRB 50, 4954 (1994)].
c
c input
c   nup    electron density [bohr**(-3)]
c   ndn    electron density [bohr**(-3)]
c   gup    abs of gradient of upspin density, |grad(nup)| [bohr**(-4)]
c   gdn    abs of gradient of downspin density, |grad(ndn)| [bohr**(-4)]
c
c output
c   fxc       exchange-correlation energy density [hartree]
c   dfxcup    d(fxc)/dnup [hartree]
c   dfxcdn    d(fxc)/dndn [hartree]
c   dfxcdgup  d(fxc)/d(|grad(nup)|) [hartree * bohr] 
c   dfxcdgdn  d(fxc)/d(|grad(ndn)|) [hartree * bohr] 
c
c Note that dfxcdgtot = d(fxc)/d(|grad(nup+ndn)|) = 0      
c
c**********************************************************************

      subroutine am05wbs(nup, ndn, gup, gdn, fx, fc, 
     .  dfxup, dfxdn, dfcup, dfcdn, dfxdgup, dfxdgdn,
     .  dfcdgup, dfcdgdn)

      implicit none

c     ** Input parameters
      real*8 nup, gup
      real*8 ndn, gdn

c     ** Output parameters
      real*8 fx, fc, dfxup, dfxdn, dfcup, dfcdn
      real*8 dfxdgup,dfxdgdn,dfcdgup,dfcdgdn

c     ** Internal parameters
      real*8 kFup, sup
      real*8 kFdn, sdn
      real*8 ex,ec 

      integer pot
      real*8 pi
      parameter (pi = 3.141592653589793238462643383279502884197d0)

c     ** Dummy parameters (not used)
      real*8 uup, udn, tup, tdn, vxup, vxdn, vcup, vcdn, stot

      kFup = (3.0d0*pi**2*2.0d0*nup)**(1.0d0/3.0d0)
      kFdn = (3.0d0*pi**2*2.0d0*ndn)**(1.0d0/3.0d0)
      call saferecp(sup,gup,(2.0d0*kFup*nup))
      call saferecp(sdn,gdn,(2.0d0*kFdn*ndn))

c     ** Not needed input for pot=1
      pot = 1 
      stot = 0.0d0
      uup = 0.0d0
      udn = 0.0d0
      tup = 0.0d0
      tdn = 0.0d0

      call am05_xscss(nup,ndn,sup,sdn,stot,uup,udn,tup,tdn,
     *     fx,fc,vxup,vxdn,vcup,vcdn,
     *     dfxup,dfxdgup,dfcup,dfcdgup,
     *     dfxdn,dfxdgdn,dfcdn,dfcdgdn,pot)

c     ** Form the output 
c     For SIESTA don't multiply by gradient since we need to divide
c     by this quantity on return to get components.
c      dfxdgup = gup*dfxdgup
c      dfxdgdn = gdn*dfxdgdn
c      dfcdgup = gup*dfcdgup
c      dfcdgdn = gdn*dfcdgdn

      return
      end

c**********************************************************************
c am05
c Calculate the Armiento Mattsson AM05 exchange-correlation energy 
c functional and various functional derivatives for different
c potential schemes.
c
c Spin polarization type: exchange index scales with the separate 
c spin densities; correlation index and prefactor both scales with
c the separate spin densities. (xscss)
c
c Input:
c   nup     electron upspin density [bohr**(-3)]
c   ndn     electron downspin density [bohr**(-3)]
c
c   sup      scaled gradient of upspin-density
c   sdn      scaled gradient of downspin-density
c   stot     scaled gradient of total density
c
c   uup      scaled grad nup * grad | grad nup |
c   udn      scaled grad ndn * grad | grad ndn |
c
c   tup      scaled laplacian of upspin-density
c   tdn      scaled laplacian of downspin-density
c
c   pot      integer: 
c              2 = calculate potential in the traditional scheme
c                  (all input needed and all output well defined)
c              1 = calculate quantities for the White and Bird,
c                  and the Pople, Gill, and Johnson schemes for
c                  potentials (u and t and stot are never touched
c                  and vx and vc give undefined output)
c              0 = don't calculate potential (u and t and stot are never 
c                  touched and only ex and ec output are well defined)
c
c Output:
c
c   fx          exchange energy density [hartree], 
c                 Total exchange Ex = Integrate[fx]
c   fc          correlation energy density [hartree]
c                 Total correlation Ec = Integrate[fc]
c    
c if pot = 1 or 2:
c
c   dfxup     d(fx)/d(nup)
c   dfxdn     d(fx)/d(ndn)
c   dfcup     d(fc)/d(nup)
c   dfcdn     d(fc)/d(ndn)
c   dfxdgup   d(fx)/d(|grad(nup)|) * 1/|grad(nup)|
c   dfxdgdn   d(fx)/d(|grad(ndn)|) * 1/|grad(ndn)|
c   dfcdgup   d(fc)/d(|grad(nup)|) * 1/|grad(nup)|
c   dfcdgdn   d(fc)/d(|grad(ndn)|) * 1/|grad(ndn)|
c
c Note that dfxcdgtot = d(fxc)/d(|grad(nup+ndn)|) * 1/|grad(nup+ndn)| = 0
c
c if pot = 2:
c
c   vxup      upspin exchange potential
c   vxdn      downspin exchange potential
c   vcup      upspin correlation potential
c   vcdn      downspin correlation potential
c
c Citation request: when using this functional, please cite:
c "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005);
c  R. Armiento and A. E. Mattsson (unpublished)."
c
c (The first paper for the AM05 functional, the second for the
c spin-polarized version)
c
c**********************************************************************
      subroutine am05_xscss(nup,ndn,sup,sdn,stot,uup,udn,tup,tdn,
     *     fx,fc,vxup,vxdn,vcup,vcdn,
     *     dfxup,dfxdgup,dfcup,dfcdgup,
     *     dfxdn,dfxdgdn,dfcdn,dfcdgdn,pot)

      implicit none

c     ** Input parameters
      real*8 nup,ndn,sup,sdn,stot,uup,udn,tup,tdn
      integer pot

c     ** Output parameters
      real*8 fx,fc,vxup,vxdn,vcup,vcdn
      real*8 dfxup,dfxdgup,dfcup,dfcdgup
      real*8 dfxdn,dfxdgdn,dfcdn,dfcdgdn

c     ** Constants
      real*8 pi, g, a, c
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (g = 0.8098d0, a = 2.804d0)
      parameter (c = 0.7168d0)

c     ** Local variables
      real*8 s2, kF, nn, ntot
      real*8 n(2), s(2), t(2), u(2)
      real*8 dfx(2), dfxdg(2), dfc(2), dfcdg(2), vx(2), vc(2)
      real*8 exlda, vxlda, eclda, vclda(2), X(2), Xsos, Xsossos
      real*8 Hx, Hxsos, Hxsossos, Hc(2), Hcsos, Hcsossos
      real*8 F, Fsos, s2Fsossos 
      real*8 szsoz, mixder
      real*8 denom, denomsos, sdenomsoss
      real*8 zfac, zosn, w

      integer i

c     ** Initialization
      fx = 0.0d0
      fc = 0.0d0

      n(1) = nup
      n(2) = ndn
      s(1) = sup
      s(2) = sdn
      t(1) = tup
      t(2) = tdn
      u(1) = uup
      u(2) = udn

      ntot = DMAX1(n(1)+n(2),1e-16+n(2),n(1)+1e-16)

      do 10, i = 1,2

         dfx(i) = 0.0d0
         dfxdg(i) = 0.0d0
         dfc(i) = 0.0d0 
         dfcdg(i) = 0.0d0 
         vx(i) = 0.0d0
         vc(i) = 0.0d0

c        ** Avoid floating point exceptions
         if(s(i) .ge. 1.0d12) then
            X(i) = 0.0d0
            Hc(i) = g
         elseif(s(i) .le. 1.0d-30) then
            X(i) = 1.0d0
            Hc(i) = 1.0d0
         else
            X(i) = 1.0d0/(1.0d0 + a*s(i)**2)
            Hc(i) = X(i) + g*(1.0d0 - X(i))
         endif
 10   continue

c     *******************
c       LDA correlation
c     *******************
      call am05_xscss_ldapwc(n,eclda,vclda)

c     ** Loop over spin up and down
      do 20, i = 1,2

c        *************
c           Cutoffs
c        *************
         if(n(i) .le. 1.0d-16) then

c           ** specifically handle a small density in the LDA limit (s=0) n->0
c           ** A small s in AM05 is always indication that we are in the
c           ** interior LDA limit, regardless of what density we have
            if(s(i) .le. 1.0d-30) then
c The following line is the original code
c               fc = fc + n(i)*eclda
c The following line is in the form needed for SIESTA
               fc = fc + eclda
               vc(i) = vclda(i)*Hc(3-i) + eclda*(1-Hc(3-i))
               dfc(i) = vclda(i)*Hc(3-i) + eclda*(1-Hc(3-i))
               goto 20
            endif

c           ** otherwise, go on as usual but with a lowest cutoff density
            n(i) = 1.0d-16
         endif

c        ** Scaling density and gradient
c        ** Ex = 0.5 * (Ex[2*spin up n] + Ex[2*spin down n] ) 
         nn = 2.0d0*n(i)
         s2 = s(i)**2
         kF = (3.0d0*pi**2*2.0d0*n(i))**(1.0d0/3.0d0)
     
c        *******************
c           LDA exchange
c        *******************
         call am05_xscss_ldax(nn,exlda,vxlda)

c        *****************************
c            Exchange energy density
c        *****************************

c        ** Airy LAA refinement function
         call am05_xscss_lambertw(s(i)**(3.0d0/2.0d0)/sqrt(24.0d0),w)

c        ** am05_lambertw give back argument if it is < 1.0e-20
c        ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
c        ** zosn = normalized z/s
         if (s(i) < 1.0e-14) then
            zosn = 1.0d0
         else
            zosn = 24.0d0**(1.0d0/3.0d0)*w**(2.0d0/3.0d0)/s(i)
         end if
         zfac = s2*(zosn*27.0d0/32.0d0/pi**2)**2

c        ** denom = denominator of Airy LAA refinement function
         denom = 1.0d0 + c*s2*zosn*(1.0d0 + zfac)**(1.0d0/4.0d0)
         F = (c*s2 + 1.0d0)/denom
      
c        ** Exchange refinement function
         Hx = X(i) + (1.0d0 - X(i))*F

c        ** Exchange energy density, Ex = Integrate[fx]
         fx = fx + 0.5d0*nn*exlda*Hx

c        ********************
c            Correlation
c        ********************

c        ** Correlation energy density, Ec = Integrate[fc]
c The following line is the original code
c         fc = fc + 0.5d0*nn*eclda*Hc(i)
c The following line is in the form needed for SIESTA
         fc = fc + 0.5d0*eclda*Hc(i)

c        ** goto next spin if we are only calculating energies
         if(pot .eq. 0) goto 20

c        ***************************
c           Exchange derivatives for White and Bird and 
c           Pople, Gill, and Johnson schemes
c        ***************************

c        ** Interpolation index derivatives: 1/s dX/ds
         Xsos = -2.0d0*a*X(i)**2

c        ** Airy LAA refinement function derivatives, 1/s dF/ds 
c        ** szsoz = s*(dz/ds)/z
         szsoz = 1.0d0/(1.0d0 + w)
         
         Fsos = c/denom**2*(2.0d0 - zosn*
     *        ((1.0d0 - c*s2)*(1.0d0 + zfac)**(1.0d0/4.0d0) +
     *        (1.0d0 + c*s2)*(1.0d0 + 3.0d0/2.0d0*zfac)/
     *        (1.0d0 + zfac)**(3.0d0/4.0d0)*szsoz))

c        ** Refinement function derivatives, 1/s dHx/ds
c        ** We use that (1 - X) = a*X*s2
         Hxsos = (1.0d0 - X(i))*Fsos - (F - 1.0d0)*Xsos

         dfx(i) = vxlda*Hx - 4.0d0/3.0d0*exlda*s2*Hxsos
         dfxdg(i) = exlda*Hxsos/((2.0d0*kF)**2*n(i))

c        *****************************
c           Correlation derivatives for White and Bird and 
c           Pople, Gill, and Johnson schemes
c        *****************************

c        ** Correlation refinement function derivatives, 1/s dF/ds 
         Hcsos = Xsos*(1.0d0 - g)

c        ** Note: n(3-i) gives the density of the *other* spin, etc.
         dfc(i) = vclda(i)/ntot*(n(1)*Hc(1)+n(2)*Hc(2)) +
     *        eclda*n(3-i)/ntot*(Hc(i) - Hc(3-i)) - 
     *        4.0d0/3.0d0*eclda*s2*Hcsos

         dfcdg(i) = eclda*Hcsos/((2.0d0*kF)**2*n(i))

c        ** goto next spin if only doing W&B/Pople et. al.
         if(pot .eq. 1) goto 20

c        ***************************
c           Exchange potential in traditional scheme
c        ***************************
c
c        ** Interpolation index derivatives: 1/s d/ds(1/s dX/ds)                   
         Xsossos = 8.0d0*a**2*X(i)**3

c        ** Airy LAA refinement function derivatives s^2 1/s d/ds (1/s dF/ds) 
c        **  mixder = szsoz + s^2*(d^2z/ds^2)/z
         mixder = (2.0d0 - w)/(2.0d0*(1.0d0 + w)**3)

c        ** denomsos = 1/s d(denom)/ds,  sdenomsoss = s*d/ds(1/s d(denom)/ds))
         denomsos = c*zosn/(1.0d0 + zfac)**(3.0d0/4.0d0)*
     *        (1.0d0 + zfac + (1.0d0 + 3.0d0/2.0d0*zfac)*szsoz)

         sdenomsoss = c*zosn/(1.0d0 + zfac)**(7.0d0/4.0d0)*
     *        (-1.0d0 - zfac*(2.0d0 + zfac)
     *        + (1.0d0 + zfac/2.0d0*(5.0d0 + 3.0d0*zfac))*mixder
     *        + 3.0d0/2.0d0*zfac*(1.0d0 + zfac/2.0d0)*szsoz**2)

         s2Fsossos = (-4.0d0*c*s2*denom*denomsos + (c*s2 + 1.0d0)*
     *        (2.0d0*s2*denomsos**2 - denom*sdenomsoss))/
     *        denom**3

c       ** Refinement function derivatives 1/s d/ds (1/s dHx/ds) 
c       ** We use that (1 - X) = a*X*s2
         Hxsossos = - 2.0d0*Fsos*Xsos + a*X(i)*s2Fsossos -
     *        (F - 1.0d0)*Xsossos

c       ** vx formula for gradient dependent functional,
         vx(i) = vxlda*(Hx - s2*Hxsos) +
     *        exlda*((4.0d0/3.0d0*s2-t(i))*Hxsos + 
     *        (4.0d0/3.0d0*s(i)**3-u(i))*s(i)*Hxsossos)

c        *****************************
c           Correlation potential in traditional scheme
c        *****************************
         
c        ** Correlation refinement function derivatives, 1/s d/ds (1/s dHc/ds)
         Hcsossos = Xsossos*(1.0d0 - g)

c        ** vc formula for gradient dependent functional,
         vc(i) = vclda(i)*(Hc(i) - s2*Hcsos) +
     *        eclda*((4.0d0/3.0d0*s2 - t(i))*Hcsos + 
     *        (4.0d0/3.0d0*s(i)**3-u(i))*s(i)*Hcsossos) +
     *     (eclda-vclda(i))*n(3-i)/ntot*
     *        (Hc(i)-Hc(3-i)-s(i)**2*Hcsos)+
     *     (eclda-vclda(3-i))*
     *        0.5d0*
     *       (stot**2*ntot**(8.0d0/3.0d0)*0.5d0**(2.0d0/3.0d0) - 
     *          s(1)**2*n(1)**(8.0d0/3.0d0) - 
     *          s(2)**2*n(2)**(8.0d0/3.0d0))/
     *       (n(i)**(5.0d0/3.0d0))/ntot*Hcsos

 20   continue

      dfxup = dfx(1)
      dfxdgup = dfxdg(1)
      dfcup = dfc(1)
      dfcdgup = dfcdg(1)
      dfxdn = dfx(2)
      dfxdgdn = dfxdg(2)
      dfcdn = dfc(2)
      dfcdgdn = dfcdg(2)

      vxup = vx(1) 
      vxdn = vx(2)
      vcup = vc(1) 
      vcdn = vc(2)

      return
      end

c     ******************************************
c       Local density approximation exchange
c
c       input
c       n        electron density [bohr**(-3)]
c
c       output
c       ex       exchange energy per electron [hartree]
c       vx       exchange potential [hartree]
c
c       Copyright (c) 2005, Rickard Armiento
c     ******************************************

      subroutine am05_xscss_ldax(n,ex,vx)
c     ** Input parameters
      real*8 n

c     ** Output parameters
      real*8 ex, vx

c     ** Constants
      real*8 pi
      parameter (pi = 3.141592653589793238462643383279502884197d0)

      vx = -(3.0d0*n/pi)**(1.0d0/3.0d0)
      ex = (3.0d0/4.0d0*vx)

      return
      end

c     ***********************************************
c       Local density approximation correlation
c
c       input
c       n(2)     electron upspin/downspin density [bohr**(-3)]
c
c       output
c       ec       correlation energy per electron [hartree]
c       vc(2)    correlation upspin/downspin potential [hartree]
c
c       As parameterized by Perdew Wang,
c         Phys. Rev. B 45, 13244 (1992) 
c       Based on Monte Carlo data by Ceperley Alder, 
c         Phys. Rev. Lett. 45, 566 (1980)
c
c       (Clean room implementation from paper)
c
c       Copyright (c) 2005, Rickard Armiento
c     ***********************************************
      subroutine am05_xscss_ldapwc(n,ec,vc)
      implicit none
c     ** Input parameters
      real*8 n(2)

c     ** Output parameters
      real*8 ec, vc(2)

c     ** Constants
      real*8 pi
      real*8 A0,a01,b01,b02,b03,b04
      real*8 A1,a11,b11,b12,b13,b14
      real*8 Aa,aa1,ba1,ba2,ba3,ba4
      parameter (pi = 3.141592653589793238462643383279502884197d0)
      parameter (a01 = 0.21370d0)
      parameter (b01 = 7.5957d0)
      parameter (b02 = 3.5876d0)
      parameter (b03 = 1.6382d0)
      parameter (b04 = 0.49294d0)
      parameter (a11 = 0.20548d0)
      parameter (b11 = 14.1189d0)
      parameter (b12 = 6.1977d0)
      parameter (b13 = 3.3662d0)
      parameter (b14 = 0.62517d0)
      parameter (aa1 = 0.11125d0)
      parameter (ba1 = 10.357d0)
      parameter (ba2 = 3.6231d0)
      parameter (ba3 = 0.88026d0)
      parameter (ba4 = 0.49671d0)
c     ** Paper actually use this:
c      parameter (A0 = 0.031091d0)
c      parameter (A1 = 0.015545d0)
c      parameter (Aa = 0.016887d0)
c     ** But routines now "defacto standard" was distributed using:
      parameter (A0 = 0.0310907d0)
      parameter (A1 = 0.01554535d0)
      parameter (Aa = 0.0168869d0)

c     ** Local variables
      real*8 xi, rsq, f, fp, fb0, mac, ec0, ec1, ecrs, ecxi
      real*8 Q0, Q1, Q1p, ec0rs, ec1rs, acrs, fdenom

c     ** Actual values from paper
c     fdenom = (2.0d0**(4.0d0/3.0d0)-2.0d0)
c     fb0 = 4.0d0/(9.0d0*(2.0d0**(1.0d0/3.0d0)-1.0d0))
c     ** Replaced with "defacto standard" approximations
      fdenom = 0.5198421d0
      fb0 = 1.709921d0

c     ** Cutoff
      if((n(1)+n(2)) .le. 1e-30) then
         ec = 0d0
         vc(1) = 0d0
         vc(2) = 0d0
         return
      endif

      xi = (n(1)-n(2))/(n(1)+n(2))
      rsq = (3.0d0/(4.0d0*pi*(n(1)+n(2))))**(1.0d0/6.0d0)
      f = ((1.0d0+xi)**(4.0d0/3.0d0)+(1.0d0-xi)**(4.0d0/3.0d0)-2.0d0)/
     *     fdenom

      mac = -2.0d0*Aa*(1.0d0 + aa1*rsq**2)* 
     *     log(1.0d0 + 1.0d0/
     *     (2.0d0*Aa*rsq*(ba1 + rsq*(ba2 + rsq*(ba3 + ba4*rsq)))))

      ec0 = -2.0d0*A0*(1.0d0 + a01*rsq**2)* 
     *     log(1.0d0 + 1.0d0/
     *     (2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))))

      ec1 = -2.0d0*A1*(1.0d0 + a11*rsq**2)* 
     *     log(1.0d0 + 1.0d0/
     *     (2.0d0*A1*rsq*(b11 + rsq*(b12 + rsq*(b13 + b14*rsq)))))

      ec = ec0 - mac*f/fb0*(1.0d0-xi**4) + (ec1-ec0)*f*xi**4

      Q0 = -2.0d0*A0*(1.0d0 + a01*rsq**2)
      Q1 = 2.0d0*A0*rsq*(b01 + rsq*(b02 + rsq*(b03 + b04*rsq)))
      Q1p = A0*(b01/rsq+2.0d0*b02+3.0d0*b03*rsq+4.0d0*b04*rsq**2)
      ec0rs = -2.0d0*A0*a01*log(1.0d0 + 1.0d0/Q1)-Q0*Q1p/(Q1**2+Q1)

      Q0 = -2.0d0*A1*(1.0d0 + a11*rsq**2)
      Q1 = 2.0d0*A1*rsq*(b11 + rsq*(b12 + rsq*(b13 + b14*rsq)))
      Q1p = A1*(b11/rsq+2.0d0*b12+3.0d0*b13*rsq+4.0d0*b14*rsq**2)
      ec1rs = -2.0d0*A1*a11*log(1.0d0 + 1.0d0/Q1)-Q0*Q1p/(Q1**2+Q1)

      Q0 = -2.0d0*Aa*(1.0d0 + aa1*rsq**2)
      Q1 = 2.0d0*Aa*rsq*(ba1 + rsq*(ba2 + rsq*(ba3 + ba4*rsq)))
      Q1p = Aa*(ba1/rsq+2.0d0*ba2+3.0d0*ba3*rsq+4.0d0*ba4*rsq**2)
      acrs = -2.0d0*Aa*aa1*log(1.0d0 + 1.0d0/Q1)-Q0*Q1p/(Q1**2+Q1)

      ecrs = ec0rs*(1.0d0 - f*xi**4) + 
     *     ec1rs*f*xi**4 - acrs*f/fb0*(1.0d0 - xi**4)
      fp = 4.0d0/3.0d0*((1.0d0+xi)**(1.0d0/3.0d0) - 
     *     (1.0d0-xi)**(1.0d0/3.0d0))/fdenom
      ecxi = 4.0d0*xi**3*f*(ec1-ec0+mac/fb0)+fp*(xi**4*ec1-xi**4*ec0+
     *     (1.0d0-xi**4)*(-mac)/fb0)

      vc(1)=ec - rsq**2/3.0d0*ecrs - (xi-1.0d0)*ecxi
      vc(2)=ec - rsq**2/3.0d0*ecrs - (xi+1.0d0)*ecxi
      
      end


c     ***********************************************
c       LambertW function. 
c
c       Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), 
c         Adv. in Comp. Math. 5(4):329-359. 
c       Implementation approach loosely inspired by the 
c       GNU Octave version by N. N. Schraudolph, but this 
c       implementation is only for real values and 
c       principal branch.
c
c       Copyright (c) 2005, Rickard Armiento
c     ***********************************************

      subroutine am05_xscss_lambertw(z,result)
      implicit none

c     input
      real*8 z
c     output
      real*8 result
c     local variables
      real*8 e,t,p    
      integer i

c     ** If z too low, go with the first term of the power expansion, z
      if( z .lt. 1.0d-20) then
         result = z
         return
      endif

      e = exp(1.0d0)

c     ** Inital guess
      if( abs(z + 1.0d0/e) .gt. 1.45d0 ) then
c        ** Asymptotic expansion at 0 and Inf
         result = log(z)
         result = result - log(result)
      else
c        ** Series expansion about -1/e to first order
         result = 1.0d0*sqrt(2.0d0*e*z + 2.0d0) - 1.0d0
      endif

c     ** Find result through iteration
      do i=1,10
         p = exp(result)
         t = result*p - z
         if( result .ne. -1.0d0 ) then
            t = t/(p*(result + 1.0d0) - 
     *           0.5d0*(result + 2.0d0)*t/(result + 1.0d0))
         else
            t = 0.0d0
         endif
         result = result - t
         if(abs(t) < (2.48d0*1.0d-14)*(1.0d0 + abs(result))) then
            return
         endif
      enddo

c     ** This should never happen!;
      write(0,*) 'am05_xscss_lambertw: iteration limit reached.' 
      write(0,*) 'Likely caused by invalid input, execution aborted.'
      stop
      return
      end

