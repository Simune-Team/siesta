c
c $Id: gexcorr.f,v 1.2 1997/05/22 17:32:14 wdpgaara Exp $
c
c $Log: gexcorr.f,v $
c Revision 1.2  1997/05/22 17:32:14  wdpgaara
c Moving from RCSfiles to ATM_1_0
c
c Revision 1.1.1.1  1997/01/07 08:38:54  wdpgaara
c Froyen-Troullier-Martins-AG atom code
c
c Revision 1.1  1991/12/14  00:34:49  alberto
c Initial revision
c
      subroutine gexcorr(id,cdd,cdu,cdc,vod,vou,vxc,vc,gga_exch_corr,
     &                   gga_correlation)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
c
c     Compute the GGA exchange and correlation energy and potential.
c     Note: This subroutine is not suitable for spin-polarized
c     calculations. It should be easy to extend to handle that case.
c     The relativistic correction to exchange has also been removed,
c     until somebody comes up with a background reference.
c
c     Revised by Alberto Garcia
c
cJLM The only major modification is that the constants for the
c    ceperly-alder 'ca' method are placed in parameter
c    statements, this was done so non-opt compiliers
c    would minimize the number of calculations.
c
C     .. Parameters ..
c
      double precision tiny_charge
      parameter (tiny_charge=1.d-12)
c
      double precision zero, one, pfive, opf
      parameter (zero=0.D0,one=1.D0,pfive=.5D0,opf=1.5D0)
      double precision psevf, c0504
      parameter (psevf=0.75D0,c0504=0.0504D0)
      double precision c0254, c014, c0406
      parameter (c0254=0.0254D0,c014=0.014D0,c0406=0.0406D0)
      double precision c15p9, c0666, c11p4
      parameter (c15p9=15.9D0,c0666=0.0666D0,c11p4=11.4D0)
      double precision c045, c7p8, c88, c20p592
      parameter (c045=0.045D0,c7p8=7.8D0,c88=0.88D0,c20p592=20.592D0)
      double precision c3p52, c0311, c0014
      parameter (c3p52=3.52D0,c0311=0.0311D0,c0014=0.0014D0)
      double precision c0538, c0096, c096
      parameter (c0538=0.0538D0,c0096=0.0096D0,c096=0.096D0)
      double precision c0622, c004, c0232
      parameter (c0622=0.0622D0,c004=0.004D0,c0232=0.0232D0)
      double precision c1686, c1p3981, c2611
      parameter (c1686=0.1686D0,c1p3981=1.3981D0,c2611=0.2611D0)
      double precision c2846, c1p0529, c3334
      parameter (c2846=0.2846D0,c1p0529=1.0529D0,c3334=0.3334D0)
      double precision con1, con2, con3
      parameter (con1=1.D0/6,con2=0.008D0/3,con3=0.3502D0/3)
      double precision con4, con5, con6
      parameter (con4=0.0504D0/3,con5=0.0028D0/3,con6=0.1925D0/3)
      double precision con7, con8, con9
      parameter (con7=0.0206D0/3,con8=9.7867D0/6,con9=1.0444D0/3)
      double precision con10, con11
      parameter (con10=7.3703D0/6,con11=1.3336D0/3)
C     ..
C     .. Scalar Arguments ..
      double precision vxc, vc, gga_exch_corr, gga_correlation
      character id*1
C     ..
C     .. Array Arguments ..
      double precision cdd(nrmax), cdu(nrmax), cdc(nrmax), vod(nrmax),
     &                 vou(nrmax)
c
      double precision dens(nrmax), gradn(nrmax), dum(nrmax),
     &                 wspl(3*nrmax), lda_xpot(nrmax), lda_cpot(nrmax),
     &                 gga_d_xen_dn(nrmax), gr_d_cen_dn(nrmax),
     &                 fun_x(nrmax), fun_c(nrmax), div_x(nrmax),
     &                 div_c(nrmax)
c
C     .. Local Scalars ..
      double precision a0, alp, be, cdsum, ftrd, pi, rs, sph_factor,
     &                 spw, rslog, sb, sqrs, te, third, s0, ax, phi,
     &                 strd, strd2, dcndn, dphidn, ro, lda_corr,
     &                 lda_exchange, gga_exchange, lda_correlation,
     &                 lda_xener_dens, lda_cener_dens,
     &                 gga_xener_dens, gr_cener_dens, gga_cener_dens,
     &                 gga_xpot, gga_cpot, vxc_adj, vc_adj
      integer i, ierr, ll
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, atan, log, sqrt
C     ..
      logical leqi
      external leqi
C     ..
      include 'gga.h'
c
      pi = 4*atan(one)
      third = one/3
c
      s0 = 1.D0/(2.D0*(3*pi*pi)**third)
      ax = -3.d0/2.d0*(3.d0/pi)**third
c
      ftrd = 4*third
      a0 = (4/(9*pi))**third
c
      alp = 2*third
c
      strd = 7.d0/3.d0
      strd2 = strd/2.d0
c
      if (leqi(id,'s')) then
c
         write(6,'(/,1x,a,/)')
     &     'GGA exch-corr not available for the spin-polarized case'
c
         stop
c
      end if
c
c     Compute the gradient of n (spatial charge dens)
c
      do 10 i = 2, nr
c
         cdsum = cdd(i) + cdu(i)
         if (ifcore .ge. 1) cdsum = cdsum + cdc(i)
c
c        Remember that cd{u,d} contain 4pir^2 rho(r)...
c
         dens(i) = cdsum/(4*pi*r(i)**2)
c
   10 continue
      dens(1) = dens(2) - (dens(3)-dens(2))*r(2)/(r(3)-r(2))
c
      call splift(r,dens,gradn,dum,nr,wspl,ierr,0,zero,zero,zero,zero)
c
c     Initialize accumulators
c
      vxc = zero
      vc = zero
      lda_exchange = zero
      lda_correlation = zero
      gga_exchange = zero
      gga_correlation = zero
c
c      start loop (at the second point...see below)
c
      ll = 4
      do 20 i = 2, nr
c
         sph_factor = 4*pi*r(i)**2
c
         fun_x(i) = 0.d0
         fun_c(i) = 0.d0
         lda_cpot(i) = 0.d0
         gga_xener_dens = 0.d0
         gga_cener_dens = 0.d0
         gga_d_xen_dn(i) = 0.d0
         gr_d_cen_dn(i) = 0.d0
c
         ro = dens(i)
c
cag  Careful with this !!
cag
cagcag         if (ro .le. tiny_charge) go to 100
c
         rs = (3.d0/(4.d0*pi*ro))**third
c
c      exchange
c
         lda_xpot(i) = -3*alp/(pi*a0*rs)
         lda_xener_dens = 3.d0*ro*lda_xpot(i)/4.d0
c
c           Exchange energy in the GGA
c
         spw = s0*abs(gradn(i))/ro**ftrd
         gga_xener_dens = f(spw)*lda_xener_dens
c
c
c        Correlation (Ceperly-Alder only)
c
c          The Perdew-Zunger parameterization is used.
c          See Phys. Rev. B 23 5075 (1981).
c
         if (rs .gt. one) then
            sqrs = sqrt(rs)
            te = one + con10*sqrs + con11*rs
            be = one + c1p0529*sqrs + c3334*rs
            lda_corr = -c2846/be
            lda_cpot(i) = lda_corr*te/be
         else
            rslog = log(rs)
            lda_corr = (c0622+c004*rs)*rslog - c096 - c0232*rs
            lda_cpot(i) = (c0622+con2*rs)*rslog - con3 - con4*rs
         end if
c
         lda_cener_dens = ro*lda_corr
c
c           Gradient contribution to the correlation energy
c
         phi = b_phi/c(rs)*abs(gradn(i))/ro**strd2
         gr_cener_dens = abs(gradn(i))**2*exp(-phi)*c(rs)/ro**ftrd
c
c           Add it to the LDA part to get the GGA
c
         gga_cener_dens = lda_cener_dens + gr_cener_dens
c
c        The exchange-correlation potential in the GGA is
c        made of two parts:
c
c        Vxc = dexc/dn - div ( 1/|delta n| dexc/d |delta n| * delta n)
c
c        Now we compute dex/dn (gga_d_xen_dn)
c        and dec/dn   (gr_d_cen_dn)    (See notes for details)
c
c        exchange:
c
         gga_d_xen_dn(i) = ftrd*ax*ro**third*f(spw)*(1.d0-g(spw))
c
c        correlation
c
         dcndn = dc(rs)*(-rs/3.D0/ro)
         dphidn = -phi*(dcndn/c(rs)+strd2/ro)
c
         gr_d_cen_dn(i) = gr_cener_dens*(-dphidn+dcndn/c(rs)-ftrd/ro)
c
c        Now calculate "divergence" terms associated with the dependence
c        on |delta n| .
c
c        Evaluate the functions whose "divergence" we need
c
c        exch : 1/|delta n| dex/d |delta n| * (delta n)
c        corr : 1/|delta n| dec/d |delta n| * (delta n)
c                                                               
c        (Careful, we have to take into account the sign of delta n)
c
c      exchange
c
         fun_x(i) = gradn(i)*gga_xener_dens*s_over_deln(ro)**2*
     &              h(spw) * r(i)**2
c
c          correlation  (gradient part only)
c
         fun_c(i) = gradn(i)*(2.D0-phi)*exp(-phi)*c(rs)/
     &              ro**ftrd * r(i)**2
c
 100     continue
c
c        Add to the integrated values. Careful with 4pir^2 factors...
c
         gga_exchange = gga_exchange + ll*sph_factor*gga_xener_dens*
     &                  rab(i)
         lda_exchange = lda_exchange + ll*sph_factor*lda_xener_dens*
     &                  rab(i)
         gga_correlation = gga_correlation +
     &                     ll*sph_factor*gga_cener_dens*rab(i)
         lda_correlation = lda_correlation +
     &                     ll*sph_factor*lda_cener_dens*rab(i)
         ll = 6 - ll
c
   20 continue
c
c     Complete Simpson's rule by multiplying by 1/3 (note that "h"
c     is 1 since we are integrating over the "i" grid)
c
      gga_exchange = gga_exchange/3
      lda_exchange = lda_exchange/3
      gga_correlation = gga_correlation/3
      lda_correlation = lda_correlation/3
c
      gga_exch_corr = gga_exchange + gga_correlation
c
c
c     Now for the potentials.
c
c     Take the derivative for the "divergence" terms
c
      fun_x(1) = fun_x(2) - (fun_x(3)-fun_x(2))*r(2)/(r(3)-r(2))
      fun_c(1) = fun_c(2) - (fun_c(3)-fun_c(2))*r(2)/(r(3)-r(2))
c
      call splift(r,fun_x,div_x,dum,nr,wspl,ierr,0,zero,zero,zero,
     &            zero)
      call splift(r,fun_c,div_c,dum,nr,wspl,ierr,0,zero,zero,zero,
     &            zero)
c
c     Final pass
c
      ll = 4
      do 30 i = 2, nr
c
         gga_xpot = gga_d_xen_dn(i) - div_x(i)/r(i)**2
         gga_cpot = lda_cpot(i) + gr_d_cen_dn(i) - div_c(i)/r(i)**2
c
         vxc_adj = gga_xpot + gga_cpot
         vc_adj = gga_cpot
c
         vod(i) = vod(i) + gga_xpot + gga_cpot
         vou(i) = vod(i)
c
c        Add to the integrated values.
c        Note that the core charge is not included for the correction
c        terms.
c
         vxc = vxc + ll*(cdd(i)+cdu(i))*vxc_adj*rab(i)
         vc = vc + ll*(cdd(i)+cdu(i))*vc_adj*rab(i)
c
         ll = 6 - ll
   30 continue
c
c     Complete Simpson's rule by multiplying by 1/3 (note that "h"
c     is 1 since we are integrating over the "i" grid)
c
      vxc = vxc/3
      vc = vc/3
c
c     Extrapolate backwards for the first point...
c
      vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
      vou(1) = vod(1)
c
      return
c
      end
