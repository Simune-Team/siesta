C $Id: cgwf.f,v 1.7 1999/02/26 14:37:06 wdpgaara Exp $

      subroutine cgwf(iscf,itmax,ftol,eta,enum,nbasis,nbands,
     .                nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .                g,hg,xi,
     .                fe,iter,dm,edm)
C *******************************************************************
C Conjugate gradient minimization of the Kim et al. functional
C Returns the converged WF. coefficients and the density matrices.
C (PRB 52, 1640 (95))
C Written by P.Ordejon, October'96
C ************************** INPUT **********************************
C integer iscf                : SCF Iteration cycle (used to find 
C                                control vectors only if iscf=1)
C integer itmax               : Maximum number of CG iterations
C real*8 ftol                 : Relative tolerance in CG minimization
C                                (recomended: 1e-8)
C real*8 eta                  : Fermi level parameter of Kim et al.
C real*8 enum                 : Total number of electrons
C integer nbasis              : Number of atomic orbitals
C integer nbands              : Number of Localized Wave Functions
C integer nhmax               : First dimension of listh and H, and maximum
C                               number of nonzero elements of each row of H
C integer numh(nbasis)        : Control vector of H matrix
C                               (number of nonzero elements of each row of H)
C integer listh(nhmax,nbasis) : Control vector of H matrix
C                               (list of nonzero elements of each row of H)
C integer ncmax               : First dimension of listc and C, and maximum
C                               number of nonzero elements of each row of C
C integer numc(nbasis)        : Control vector of C matrix
C                               (number of nonzero elements of each row of C)
C integer listc(ncmax,nbasis) : Control vector of C matrix
C                               (list of nonzero elements of each row of C)
C real*8 h(nhmax,nbasis)      : Hamiltonian matrix (sparse)
C real*8 s(nhmax,nbasis)      : Overlap matrix (sparse)
C real*4 g(ncmax,nbasis)      : Auxiliary space matrix
C real*4 hg(ncmax,nbasis)     : Auxiliary space matrix
C real*8 xi(ncmax,nbasis)     : Auxiliary space matrix (gradients)
C ******************** INPUT AND OUTPUT *****************************
C real*8 c(ncmax,nbasis)      : Localized Wave Functions (sparse)
C ************************* OUTPUT **********************************
C real*8 fe                   : Final electronic band structure energy
C integer iter                : Number of CG iterations
C real*8 dm(nhmax,nbasis)     : Density matrix (sparse)
C real*8 edm(nhmax,nbasis)    : Energy Density matrix (sparse)
C *******************************************************************
      implicit none

      integer
     .  nbands,nbasis,ncmax,nhmax

      integer
     .  iscf,iter,itmax,
     .  listc(ncmax,nbasis),listh(nhmax,nbasis),
     .  numc(nbasis),numh(nbasis)

      real*4
     .  g(ncmax,nbasis),hg(ncmax,nbasis)

      double precision
     .  c(ncmax,nbasis),dm(nhmax,nbasis),edm(nhmax,nbasis),
     .  enum,eta,fe,ftol,
     .  h(nhmax,nbasis),s(nhmax,nbasis),xi(ncmax,nbasis)

C  Internal variables .................................................
      integer 
     . i,its,j,iopt,irestart,numit

      double precision
     .  dgg,e,e3(3),elamb,eps,fp,gam,gg,lam,lambda,lm,partial

      logical 
     .  iout,itest

      parameter (eps=1.d-15)
      parameter (irestart=300)
C ..........................


C Set step for line minim: lam (empirically 1.0e-1 works fine) 
      lam = 1.0d0
      lm = 0.0d0
      iter = 0
      partial = 0.0d0
      itest = .false.
      numit = 0

      if (iter .eq. itmax) goto 15

C Calculate control vectors only if in first SCF step ....................
      if (iscf .eq. 1) then
        iopt = 0
        call eandg(iopt,eta,enum,lam,
     .             nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .             nbasis,nbands,
     .             e3,e,xi,dm,edm)
      endif
C .....................

C Calculate gradient for first CG step ...................................
      iopt = 2
      call eandg(iopt,eta,enum,lam,
     .           nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .           nbasis,nbands,
     .           e3,fp,xi,dm,edm)
C .....................

      elamb = fp

      do 11 i = 1,nbasis
      do 11 j = 1,numc(i)
        g(j,i) = -xi(j,i)
        hg(j,i) = g(j,i)
        partial = partial + hg(j,i) * xi(j,i)
        xi(j,i) = hg(j,i)
11    continue


C Loop for the CG minimization ____________________________________________
      do 14 its = 1,itmax
        iter = its
        numit = numit + 1

        iout = .true.
        
        write(6,33) its,partial,fp

C Check if the gradient at current point is negative, and sufficiently
C large to avoid problems. Otherwise, set a fixed lambda
C to CG iteration cycle with iout = .false.  ................................
        if (abs(partial) .le. 1.5e-8) then
          lambda = lm*0.4
          if (lm .eq. 0.0) lambda=0.4
C         iout = .false.
          goto 1000
        endif
        if (partial .gt. 0.0e0) then
          lambda = -0.0001
          iout = .false.
          goto 1000
        endif
C ...........................

C Get the energy at three points separated by lam ...........................
        iopt = 1
        call eandg(iopt,eta,enum,lam,
     .             nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .             nbasis,nbands,
     .             e3,e,xi,dm,edm)
C ...........................

C Solve exactly the line minimzation with a forth order polinomial ..........
        call linmin4(partial,elamb,e3,lam,lambda)
C ...........................

C Check if the solution of the line minimization is reasonable 
C (ie, lambda is small)
C Otherwise, set a fixed lambda and return to CG iteration 
C cycle with iout = .false.  ................................
        if (abs(lambda) .ge. 10.0e0) then
          lambda = lm * 0.3
          if (lm .eq. 0.0) lambda = 0.01
          iout = .false.
          goto 1000
        endif
C ...........................

        iout = .true.

1000    continue

C Update point to line mimimum and multiply gradient by lambda .............
        do i = 1,nbasis
          do j = 1,numc(i)
            xi(j,i) = lambda * xi(j,i)
            c(j,i) = c(j,i) + xi(j,i)
          enddo
        enddo
C ...........................

C Calculate average lambda (lm) during CG minimization .....................
        lm = float(its-1) / float(its) * lm + lambda / float(its)
C ...........................


C Calculate energy and gradient at current point ...........................
        partial=0.0d0
        iopt = 2
        call eandg(iopt,eta,enum,lam,
     .             nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .             nbasis,nbands,
     .             e3,fe,xi,dm,edm)
C ...........................

C Check if minimization has converged ......................................
        if (2.*abs(fe-fp).le.ftol*(abs(fe)+abs(fp)+eps)) then
          if (iout) then
          write(6,"(/a)") 'cgwf:  CG tolerance reached'
          goto 16
          endif
        endif
C ...........................

C Continue if tol not reached ..............................................
        fp = fe
        elamb = fe
        gg = 0.0d0
        dgg = 0.0d0
        do 12 i = 1,nbasis
        do 12 j = 1,numc(i)
          gg = gg+g(j,i)**2
C         dgg = dgg + xi(j,i)**2
          dgg = dgg + (xi(j,i) + g(j,i)) * xi(j,i)
12      continue
        if (gg .eq. 0.0d0) goto 16
        gam = dgg / gg
        do 13 i = 1,nbasis
        do 13 j = 1,numc(i)
          g(j,i) = -xi(j,i)
          if (itest) gam = 0.0d0
          hg(j,i) = g(j,i) + gam * hg(j,i)
          partial = partial + hg(j,i) * xi(j,i)
          xi(j,i) = hg(j,i)
13      continue

      itest = .false.
      if (numit .eq. irestart) then
        itest = .true.
        numit = 0
      endif

14    continue
C ......................
C end CG loop ________________________

C Exit if maximum number of iterations has been reached ...................
15    continue

      write(6,"(/a)") 'cgwf: Maximum number of CG iterations reached'

16    continue

C Compute density matrix ...................................................
          iopt = 3
          call eandg(iopt,eta,enum,lam,
     .               nhmax,numh,listh,ncmax,numc,listc,h,s,c,
     .               nbasis,nbands,
     .               e3,fe,xi,dm,edm)
C ............................
C33   format('ITER = ',i4,6x,'GRAD = ',f18.6,6x,'ENER = ',f14.6)
 33   format('cgwf: iter = ',i4,6x,'grad = ',f18.6,6x,'Eb(Ry) = ',f14.6)
      return
      end
C ......................
C *******************************************************************



      subroutine linmin4(partial,elamb,ener,lam,lambda)
C *******************************************************************
C  Subroutine for the exact minimization of a 4th order polinomial.
C  Uses the values of the polinomial at four points, and the value of
C  the derivative at one point to determine the polynomial.
C  From the polinomial coefficients, it determines the minimum.
C
C  The input are the values of the polinomial and its derivative
C  at a point x0, and the values of the polinomial at three other
C  points x_i = x0 + lam * i.
C 
C  The output is the distance lamba from the initial point x0 to the 
C  minimum x = x0 + lambda
C
C  P.Ordejon, 12/92 - 4/93. Re-written 10/96
C *********************** INPUT *************************************
C real*8 partial              : Derivative at current point
C real*8 elamb                : Value of polinomial at curent point
C real*8 ener(3)              : Value of polinomial at three points
C real*8 lam                  : Step between points
C ********************** OUTPUT *************************************
C real*8 lambda               : Distance to minimum
C *******************************************************************
      implicit none

      double precision
     .  elamb,partial,lam,lambda,ener(3)

C Internal variables ......................................................
      integer 
     .  in,i,j,jj

      double precision 
     .  b(5,5),bi(5,5),e(5),e0,e1,e2,e3,ep1,ep2,ep3,
     .  l,lambda1,lambda2,lambda3,q,r,sol1,sol2,test

      complex 
     .  c


C Set up linear equations system to calculate polinomial coeffs.............
      e0 = 0.0d0
      e1 = 0.0d0
      e2 = 0.0d0
      e3 = 0.0d0

      e(1) = elamb
      e(2) = partial
      b(1,1) = 1.0d0
      b(1,2) = 0.0d0
      b(1,3) = 0.0d0
      b(1,4) = 0.0d0
      b(1,5) = 0.0d0
      b(2,1) = 0.0d0
      b(2,2) = 1.0d0
      b(2,3) = 0.0d0
      b(2,4) = 0.0d0
      b(2,5) = 0.0d0

      do in = 1,3
        l = float(in)*lam
        B(in+2,1) = 1.0d0
        do j = 2,5
          B(in+2,j) = l**(j-1)
        enddo
      enddo


      do i = 1,3
        e(i+2) = ener(i)
      enddo

      call inver(b,bi,5)

      do jj = 1,5
        e0 = e0 + bi(2,jj) * e(jj)
        e1 = e1 + bi(3,jj) * e(jj) * 2.0
        e2 = e2 + bi(4,jj) * e(jj) * 3.0
        e3 = e3 + bi(5,jj) * e(jj) * 4.0
      enddo
C .....................


C Solve the minimum ......................................................
      q = (1.d0 / 3.d0) * (e1 / e3) - 
     .    (1.d0 / 9.d0) * (e2 / e3)**2
      r = (1.d0 / 6.d0) * ((e1 * e2) / (e3 * e3) - 
     .     3.d0 * (e0 / e3)) - (1.d0 / 27.d0) * ( e2/e3 )**3

      test = q**3 + r**2

      if (test .le. 0.0d0) goto 1234

      if (r + test**0.5d0 .lt. 0.0d0) then
        sol1 = - ((-r - test**0.5d0)**(1.d0 / 3.d0))
        goto 122
      endif
      sol1 = (r + test**0.5d0)**(1.d0 / 3.d0)
122   continue
      if (r .lt. test**0.5d0) then
        sol2 = - ( (-r + test**0.5d0)**(1.d0 / 3.d0))
        goto 123
      endif
      sol2 = (r - test**0.5d0)**(1.d0 / 3.d0)
123   continue

      lambda = (sol1 + sol2) - (e2 / e3) / 3.0d0

      goto 1235

1234  continue

      c = (0.0,1.0) * (-test)**0.5
      c = c + r * (1.0,0.0)
      c = cexp(1.e0 / 3.e0 * clog(c))

      lambda1 = 2.0d0 * real(c) - (e2 / e3) / 3.0d0

      lambda2 = (-(e2 / e3 + lambda1) + dsqrt((e2 / e3 + lambda1)**2 +
     .          4.d0 * (e0 / e3) / lambda1)) / 2.0d0

      lambda3 = (-(e2 / e3 + lambda1) - dsqrt((e2 / e3 + lambda1)**2 +
     .          4.d0 * (e0 / e3) / lambda1)) / 2.0d0

      ep1 = lambda1 * e0 + lambda1**2 * e1 / 2.0d0 + 
     .      lambda1**3 * e2 / 3.0d0 + 
     .      lambda1**4 * e3 / 4.0d0

      ep2 = lambda2 * e0 + lambda2**2 * e1 / 2.0d0 +
     .      lambda2**3 * e2 / 3.0d0 +
     .      lambda2**4 * e3 / 4.0d0

      ep3 = lambda3 * e0 + lambda3**2 * e1 / 2.0d0 + 
     .      lambda3**3 * e2 / 3.0d0 +
     .      lambda3**4 * e3 / 4.0d0

      if (ep1 .lt. ep2) then
        lambda=lambda1
        goto 4321
      endif
      lambda=lambda2
4321  continue
      if ((ep3 .lt. ep1) .and. (ep3 .lt. ep2)) then
        lambda=lambda3
      endif
C .....................

1235  continue


      return
      end

