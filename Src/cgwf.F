      subroutine cgwf(iscf,itmax,ftol,eta,enum,nbasis,nbands,
     .                maxnh,numh,listhptr,listh,ncmax,numc,
     .                listc,h,s,c,g,hg,xi,fe,iter,dm,edm,nbasisloc,
     .                maxnf,maxnft,maxnhij,maxnct)
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
C integer maxnh               : First dimension of listh and H, and maximum
C                               number of nonzero elements of H
C integer numh(nbasisloc)     : Control vector of H matrix
C                               (number of nonzero elements of each row of H)
C integer listhptr(nbasisloc) : Control vector of H matrix
C                               (pointer to each row of H matrix)
C integer listh(maxnh)        : Control vector of H matrix
C                               (list of nonzero elements of each row of H)
C integer ncmax               : First dimension of listc and C, and maximum
C                               number of nonzero elements of each row of C
C integer numc(nbasis)        : Control vector of C matrix
C                               (number of nonzero elements of each row of C)
C integer listc(ncmax,nbasis) : Control vector of C matrix
C                               (list of nonzero elements of each row of C)
C real*8 h(maxnh)             : Hamiltonian matrix (sparse)
C real*8 s(maxnh)             : Overlap matrix (sparse)
C real*4 g(ncmax,nbasisloc)   : Auxiliary space matrix
C real*4 hg(ncmax,nbasisloc)  : Auxiliary space matrix
C real*8 xi(ncmax,nbasisloc)  : Auxiliary space matrix (gradients)
C integer nbasisloc           : Number of local basis functions stored 
C integer maxnf               : )
C integer maxnft              : )Array dimensions to be passed to lower routines
C integer maxnhij             : )
C integer maxnct              : )
C ******************** INPUT AND OUTPUT *****************************
C real*8 c(ncmax,nbasisloc)   : Localized Wave Functions (sparse)
C ************************* OUTPUT **********************************
C real*8 fe                   : Final electronic band structure energy
C integer iter                : Number of CG iterations
C real*8 dm(maxnh)            : Density matrix (sparse)
C real*8 edm(maxnh)           : Energy Density matrix (sparse)
C *******************************************************************
C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi
#endif

      implicit none

      integer
     .  nbands,nbasis,ncmax,maxnh,nbasisloc,maxnf,maxnft,maxnhij,
     .  maxnct

      integer
     .  iscf,iter,itmax,
     .  listc(ncmax,nbasis),listh(*),listhptr(nbasisloc),
     .  numc(nbasis),numh(nbasisloc)

      REAL(SP)
     .  g(ncmax,nbasisloc),hg(ncmax,nbasisloc)

      real(dp)
     .  c(ncmax,nbasisloc),dm(maxnh),edm(maxnh),enum,eta,fe,ftol,
     .  h(maxnh),s(maxnh),xi(ncmax,nbasisloc)

C  Internal variables .................................................
      integer 
     . i,ii,its,j,iopt,irestart,numit,Node,Nodes

#ifdef MPI
      integer MPIerror
      real(dp) rtmp
#endif

      real(dp)
     .  dgg,e,e3(3),elamb,eps,fp,gam,gg,lam,lambda,lm,partial

      logical 
     .  iout,itest

      parameter (eps=1.e-15_dp)
      parameter (irestart=300)
C ..........................

C Get Node number 
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C Set step for line minim: lam (empirically 1.0e-1 works fine) 
      lam = 1.0_dp
      lm = 0.0_dp
      iter = 0
      partial = 0.0_dp
      itest = .false.
      numit = 0

      if (iter .eq. itmax) goto 15

C Calculate control vectors only if in first SCF step ....................
      if (iscf .eq. 1) then
	iopt = 0
        call eandg(iopt,eta,enum,lam,
     .             maxnh,numh,listhptr,listh,ncmax,numc,listc,h,s,c,
     .             nbasis,nbands,e3,e,xi,dm,edm,nbasisloc,
     .             maxnf,maxnft,maxnhij,maxnct)
      endif
C .....................

C Calculate gradient for first CG step ...................................
      iopt = 2
      call eandg(iopt,eta,enum,lam,
     .           maxnh,numh,listhptr,listh,ncmax,numc,listc,h,s,c,
     .           nbasis,nbands,e3,fp,xi,dm,edm,nbasisloc,
     .           maxnf,maxnft,maxnhij,maxnct)
C .....................

      elamb = fp

      do i = 1,nbasisloc
        call LocalToGlobalOrb(i,Node,Nodes,ii)
        do j = 1,numc(ii)
          g(j,i) = -xi(j,i)
          hg(j,i) = g(j,i)
          partial = partial + hg(j,i) * xi(j,i)
          xi(j,i) = hg(j,i)
        enddo
      enddo

#ifdef MPI
C Global reduction of Gnorm
      call MPI_AllReduce(partial,rtmp,1,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      partial = rtmp
#endif

C Loop for the CG minimization ____________________________________________
      do 14 its = 1,itmax
        iter = its
	numit = numit + 1

        iout = .true.
        
        if (Node.eq.0) write(6,33) its,partial,fp

C Check if the gradient at current point is negative, and sufficiently
C large to avoid problems. Otherwise, set a fixed lambda
C to CG iteration cycle with iout = .false.  ................................
        if (abs(partial) .le. 1.5e-8_dp) then
  	  lambda = lm*0.4_dp
	  if (lm .eq. 0.0_dp) lambda=0.4_dp
C	  iout = .false.
	  goto 1000
        endif
        if (partial .gt. 0.0_dp) then
	  lambda = -0.0001_dp
	  iout = .false.
	  goto 1000
        endif
C ...........................

C Get the energy at three points separated by lam ...........................
        iopt = 1
        call eandg(iopt,eta,enum,lam,
     .             maxnh,numh,listhptr,listh,ncmax,numc,listc,h,s,c,
     .             nbasis,nbands,e3,e,xi,dm,edm,nbasisloc,
     .             maxnf,maxnft,maxnhij,maxnct)
C ...........................

C Solve exactly the line minimzation with a forth order polinomial ..........
        call linmin4(partial,elamb,e3,lam,lambda)
C ...........................

C Check if the solution of the line minimization is reasonable 
C (ie, lambda is small)
C Otherwise, set a fixed lambda and return to CG iteration 
C cycle with iout = .false.  ................................
        if (abs(lambda) .ge. 10.0_dp) then
          lambda = lm * 0.3_dp
          if (lm .eq. 0.0_dp) lambda = 0.01_dp
          iout = .false.
	  goto 1000
        endif
C ...........................

	iout = .true.

1000    continue

C Update point to line mimimum and multiply gradient by lambda .............
        do i = 1,nbasisloc
          call LocalToGlobalOrb(i,Node,Nodes,ii)
          do j = 1,numc(ii)
            xi(j,i) = lambda * xi(j,i)
            c(j,i) = c(j,i) + xi(j,i)
          enddo
        enddo
C ...........................

C Calculate average lambda (lm) during CG minimization .....................
        lm = real(its-1,dp) / real(its,dp) * lm
     $                  + lambda / real(its,dp)
C ...........................


C Calculate energy and gradient at current point ...........................
        partial=0.0_dp
        iopt = 2
        call eandg(iopt,eta,enum,lam,
     .             maxnh,numh,listhptr,listh,ncmax,numc,listc,h,s,c,
     .             nbasis,nbands,e3,fe,xi,dm,edm,nbasisloc,
     .             maxnf,maxnft,maxnhij,maxnct)
C ...........................

C Check if minimization has converged ......................................
        if (2.0_dp*abs(fe-fp).le.ftol*(abs(fe)+abs(fp)+eps)) then
	  if (iout) then
            if (Node.eq.0) then
	      write(6,"(/a)") 'cgwf:  CG tolerance reached'
            endif
            goto 16
	  endif
        endif
C ...........................

C Continue if tol not reached ..............................................
        fp = fe
        elamb = fe
        gg = 0.0_dp
        dgg = 0.0_dp
        do i = 1,nbasisloc
          call LocalToGlobalOrb(i,Node,Nodes,ii)
	  do j = 1,numc(ii)
            gg = gg+g(j,i)**2
            dgg = dgg + (xi(j,i) + g(j,i)) * xi(j,i)
          enddo
        enddo
#ifdef MPI
        call MPI_AllReduce(gg,rtmp,1,DAT_double,MPI_sum,
     .    MPI_Comm_World,MPIerror)
        gg = rtmp
        call MPI_AllReduce(dgg,rtmp,1,DAT_double,MPI_sum,
     .    MPI_Comm_World,MPIerror)
        dgg = rtmp
#endif
        if (gg .eq. 0.0_dp) goto 16
        gam = dgg / gg
        do i = 1,nbasisloc
          call LocalToGlobalOrb(i,Node,Nodes,ii)
	  do j = 1,numc(ii)
            g(j,i) = -xi(j,i)
	    if (itest) gam = 0.0_dp
            hg(j,i) = g(j,i) + gam * hg(j,i)
            partial = partial + hg(j,i) * xi(j,i)
            xi(j,i) = hg(j,i)
          enddo
        enddo
#ifdef MPI
        call MPI_AllReduce(partial,rtmp,1,DAT_double,MPI_sum,
     .    MPI_Comm_World,MPIerror)
        partial = rtmp
#endif

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

      if (Node.eq.0) then
        write(6,"(/a)") 'cgwf: Maximum number of CG iterations reached'
      endif

16    continue

C Compute density matrix ...................................................
      iopt = 3
      call eandg(iopt,eta,enum,lam,
     .           maxnh,numh,listhptr,listh,ncmax,numc,listc,h,s,c,
     .           nbasis,nbands,e3,fe,xi,dm,edm,nbasisloc,
     .           maxnf,maxnft,maxnhij,maxnct)
C ............................
 33   format('cgwf: iter = ',i4,6x,'grad = ',f18.6,6x,'Eb(Ry) = ',f14.6)

      end


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

      use precision

      implicit none

      real(dp), intent(in)  ::  elamb, partial, lam ,ener(3)
      real(dp), intent(out) ::  lambda

C Internal variables ......................................................

      real(dp), parameter   ::  one_third = 1.0_dp/3.0_dp

      integer in,i,j,jj

      real(dp) 
     .  b(5,5),bi(5,5),e(5),e0,e1,e2,e3,ep1,ep2,ep3,
     .  l,lambda1,lambda2,lambda3,q,r,sol1,sol2,test

      complex(dp)  c


C Set up linear equations system to calculate polinomial coeffs.............
      e0 = 0.0_dp
      e1 = 0.0_dp
      e2 = 0.0_dp
      e3 = 0.0_dp

      e(1) = elamb
      e(2) = partial
      b(1,1) = 1.0_dp
      b(1,2) = 0.0_dp
      b(1,3) = 0.0_dp
      b(1,4) = 0.0_dp
      b(1,5) = 0.0_dp
      b(2,1) = 0.0_dp
      b(2,2) = 1.0_dp
      b(2,3) = 0.0_dp
      b(2,4) = 0.0_dp
      b(2,5) = 0.0_dp

      do in = 1,3
        l = real(in,dp)*lam
        B(in+2,1) = 1.0_dp
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
        e1 = e1 + bi(3,jj) * e(jj) * 2.0_dp
        e2 = e2 + bi(4,jj) * e(jj) * 3.0_dp
        e3 = e3 + bi(5,jj) * e(jj) * 4.0_dp
      enddo
C .....................


C Solve the minimum ......................................................
      q = one_third * (e1 / e3) - 
     .    (1._dp / 9._dp) * (e2 / e3)**2
      r = (1._dp / 6._dp) * ((e1 * e2) / (e3 * e3) - 
     .     3._dp * (e0 / e3)) - (1._dp / 27._dp) * ( e2/e3 )**3

      test = q**3 + r**2

      if (test .gt. 0.0_dp) then

         if (r + sqrt(test) .lt. 0.0_dp) then
            sol1 = - ((-r - sqrt(test))**one_third)
         else
            sol1 = (r + sqrt(test))**one_third
         endif

         if (r .lt. sqrt(test)) then
            sol2 = - ( (-r + sqrt(test))**one_third)
         else
            sol2 = (r - sqrt(test))**one_third
         endif

         lambda = (sol1 + sol2) - one_third * (e2 / e3) 

      else

         c = (0.0_dp,1.0_dp) * sqrt(-test)
         c = c + r
         c = exp(one_third * log(c))

         lambda1 = 2.0_dp * real(c) - one_third * (e2 / e3)
         
         lambda2 = (-(e2/e3 + lambda1) + dsqrt((e2/e3 + lambda1)**2 +
     .        4._dp * (e0 / e3) / lambda1)) / 2.0_dp
         
         lambda3 = (-(e2/e3 + lambda1) - dsqrt((e2/e3 + lambda1)**2 +
     .        4._dp * (e0 / e3) / lambda1)) / 2.0_dp

         ep1 = lambda1 * e0 + lambda1**2 * e1 / 2.0_dp + 
     .        lambda1**3 * e2 / 3.0_dp + 
     .        lambda1**4 * e3 / 4.0_dp

         ep2 = lambda2 * e0 + lambda2**2 * e1 / 2.0_dp +
     .        lambda2**3 * e2 / 3.0_dp +
     .        lambda2**4 * e3 / 4.0_dp

         ep3 = lambda3 * e0 + lambda3**2 * e1 / 2.0_dp + 
     .        lambda3**3 * e2 / 3.0_dp +
     .        lambda3**4 * e3 / 4.0_dp

         if (ep1 .lt. ep2) then
            lambda=lambda1
         else
            lambda=lambda2
         endif

         if ((ep3 .lt. ep1) .and. (ep3 .lt. ep2)) then
            lambda=lambda3
         endif

      endif

      end

