c---
C     Parameters and Statement Functions for the 
c     Generalized Gradient Approximation
c
c     Alberto Garcia, November 1991
c
      double precision alphagr, betagr, gammagr, deltagr, lambdagr, c0,
     &                 c1, c_infinity
c
      parameter (alphagr=0.023266D0)
      parameter (betagr=7.389D-6)
      parameter (gammagr=8.723D0)
      parameter (deltagr=0.472D0)
      parameter (lambdagr=1.d4*betagr)
c
      parameter (c0=0.001667d0)
      parameter (c1=0.002568D0)
      parameter (c_infinity=2*(c0+c1))
c
c     Fgr_ca was chosen by Perdew+Wang so that the GGA correlation energy
c     would be exact for the Ne atom while using the Ceperley-Alder
c     form for the LDA correlation. Accordingly, it only works in tandem
c     with CA. Other values of fgr could in principle be calculated for
c     other LDA correlation schemes.
c
      double precision fgr_ca, b_phi
      parameter (fgr_ca = 0.11D0)
      parameter (b_phi = 1.745D0*fgr_ca*c_infinity)

c
      double precision a_f, b_f, c_f, m_f
c
      parameter (a_f=1.296d0)
      parameter (b_f=14.d0)
      parameter (c_f=0.2d0)
      parameter (m_f=1.d0/15.d0)
c
      double precision r_s, s, density
      double precision c_num, c_den, dc_num, dc_den, c, dc
      double precision beta, f, g, h
      double precision s_over_deln
c
      integer raw_int, base_int
      integer reduce
c
c     This is C(n) in J. P. Perdew: PRB 33, 8822 (1986)  eq. 6
c     BUT note that we use rydberg units instead of hartree.
c
      c_num(r_s) = c1 + alphagr*r_s + betagr*r_s**2
      dc_num(r_s) = alphagr + 2*betagr*r_s
      c_den(r_s) = (1.D0+gammagr*r_s+deltagr*r_s**2+lambdagr*r_s**3)
      dc_den(r_s) = gammagr + 2*deltagr*r_s + 3*lambdagr*r_s**2
c
      c(r_s) = 2 * ( c0 + c_num(r_s)/c_den(r_s) )
c
c     Derivative with respect to r_s
c
      dc(r_s) = 2*(dc_num(r_s)-dc_den(r_s)*c_num(r_s)/c_den(r_s))/
     &          c_den(r_s)
c
c     This is F(s) in Perdew-Wang.     G(s) = s/F * dF/ds
c
c     Introduce also H(s) such that G(s) = s**2 H(s), with H(s)
c     finite as s-->0
c
      beta(s) = 1.d0 + a_f*s**2 + b_f*s**4 + c_f*s**6
      F(s) = beta(s)**m_f
      H(s) = 2*m_f * (a_f+2*b_f*s*s+3*c_f*s**4)/beta(s)
      G(s) = s**2 * H(s)
c
      s_over_deln(density) =  density**(-4.d0/3.d0) /
     &                       (2.D0*(3*pi*pi)**(1.d0/3.d0)) 
c---
