C $Id: chempot.f,v 1.3 1999/03/11 19:26:36 ordejon Exp $

      subroutine chempot(h,s,listh,numh,rcoor,pmax,beta,
     .                   lasto,cell,xa,enum,nbasis,natoms,nhmax,
     .                   chpot,emax,emin)

C     .,gap,homo,lumo)
C *****************************************************************
C Calculates the maximum and minimum eigenvalues, the chemical 
C potential and the HOMO-LUMO gap.
C The calculation of the max. and min. eigenvalues is done with
C the Lanczos method. The chemical potential is calculated 
C with the projection method of Goedecker, and the HOMO-LUMO
C gap is obtained by the Folded-Spectrum method (combined with
C Lanczos)
C
C NOTE: In this version, the calculation of the HOMO, LUMO
C       and gap is disabled because it does work properly.
C
C Written by Maider Machado and P.Ordejon, June'98 
C ***************************** INPUT ***************************** 
C real*8 h(nhmax,nbasis)          : Hamiltonian in sparse form
C real*8 s(nhmax,nbasis)          : Overlap in sparse form
C integer listh(nhmax,nbasis)     : Control vector of sparse hamilt.
C integer numh(nbasis)            : Control vector of sparse hamilt.
C real*8 rcoor                    : Cutoff range for the projection 
C                                   vectors
C integer pmax                    : Maximum order of Chebishev expansion
C real*8 beta                     : Inverse Temperature for Chebi expansion
C integer lasto(0:natoms)         : Index vector of last orbital of 
C                                   each atom
C real*8 cell(3,3)                : Lattice vectors
C real*8 xa(3,natoms)             : Atomic positions
C real*8 enum                     : Total number of electrons
C integer nbasis                  : Number of basis orbitals
C integer natoms                  : Number of atoms
C integer nhmax                   : Maximum number of non-zero elements
C                                   of each column of sparse Hamilt.
C **************************** OUTPUT *******************************
C real*8 chpot                    : Chemical potential
C real*8 emax                     : Maximum eigenvalue
C real*8 emin                     : Minimum eigenvalue
C *************************** DISABLED ******************************
C real*8 gap                      : Energy gap
C real*8 homo                     : Highest occ. molecular orbital energy
C real*8 lumo                     : Lowest unocc. molecular orbital energy
C *********************************************************************

      implicit none

      include 'ordern.h'

      integer 
     .  natoms, nbasis, nhmax, pmax

      integer
     .  lasto(0:natoms), listh(nhmax,nbasis), numh(nbasis)

      real*8
     .  beta, cell(3,3), chpot, emin, emax, enum, 
     .  h(nhmax,nbasis),
     .  rcoor, s(nhmax,nbasis), xa(3,natoms)
C    .  homo, lumo, gap


C Internal variables...

      integer
     . ipmax, nnmax, nvmax
C nnmax = maximum number of neighbors atoms within rcoor
      parameter (nnmax = 1000)

C nvmax = maximum number of non-zero elements within a sparse vector v
      parameter(nvmax = 1000)

C ipmax = maximum order the Chebyshev polynomial exp. of the Fermi operator
      parameter(ipmax = 200)

      real*8 
     .  mu1,mu2
      parameter(mu1=-1.2,mu2=1.2)

      integer
     .  i, ia, ii, imu, iorb, indexloc(nnmax), j, ja, jan(nnmax), 
     .  ji, jj, jorb, k, listhp(maxnh,maxo), listhpp(maxnh,maxo), 
     .  listvt(maxo), listv(nvmax), m, mu, nb, nna, norb, nu, num,
     .  numhp(maxo), numloc, numv, p

      real*8
     .  betap, c(ipmax), chpotsh, deltae, emean,
     .  hbar(maxnh,maxo), Hdense(maxnh,maxnh),
     .  numb, paux(maxnh), qtot, r2ij(nnmax), 
     .  rmax, rrmod, ri(3), rr(maxo,0:ipmax), Sdense(maxnh,maxnh), tol,
     .  v(nvmax,3), vec(maxnh), xij(3,nnmax), zbrent
C    .  delta, eref2, eref1, eig2, eig1

C common to pass information to function numb ...
      common /vect/rr,c,qtot,p,nb
C ...
C common to pass beta to chebfd ...
      common /beta/betap
C ...
      external 
     .  numb, zbrent
C ...

C Assign information for Chebyshev expansion .................
C Total charge
      qtot = enum
      nb = nbasis
      p = pmax

C *** Calculate H' = S(-1)*H  (See Gibson et al, PRB 47, 9229 (92)) ***
C This is done by Choleski decomposition of S, and solving the
C linear system S H' = H, in the subspace of orbitals which overlap
C with those of a given atom. The advantage is that the Choleski
C decomposition can be done only once for each atom, and use
C the result for all the orbitas.
C Since different orbitals in the same atom have differen
C cutoff radii, this must be done carefully, since the list of
C neighbors is not the same for all orbitals.


c  Loop over atoms.............
      do ia=1,natoms

C   form S and S-1 matrices in reduced space...

C   first determine which is the longer range orbital of atom ia
C   (which will determine the reduce space) -
        num=0
        mu=0
        do i=lasto(ia-1)+1, lasto(ia)
          if(numh(i) .gt. num) then
            mu=i
            num=numh(i)
          endif
        enddo
        if (mu .eq. 0) then
          write(6,*) 'chempot: ERROR: zero neighbors for orbital ',i
          stop
        endif
        nu=numh(mu) 
C -

C Initialize overlap and interaction in reduced space -
        do i=1,nu
          do j=1,nu
            Sdense(i,j)=0.
            Hdense(i,j)=0.
          enddo
        enddo
C -

C Construct S and H in reduced space  -
C  loop over orbitals ii in reduced space
        do i=1,nu
          ii = listh(i,mu)
C  loop over orbitals jj in reduced space
          do j=i,nu
            jj = listh(j,mu)
C  see if orbitals ii and jj interact
            do k = 1,numh(ii)
              if (listh(k,ii) .eq. jj) then
                Sdense(i,j)=s(k,ii)
                Sdense(j,i)=Sdense(i,j)
                Hdense(i,j)=h(k,ii)
                Hdense(j,i)=Hdense(i,j)
              endif
            enddo
          enddo
        enddo

C  Cholesky factorization:
        call choldc(Sdense,nu,maxnh,paux)
C ...

C loop over orbitals of atom ia ....
C 
        do i=1,nu
          ii=listh(i,mu)
C  check if ii is in atom ia
          if (ii .ge. (lasto(ia-1)+1) .and. ii .le. lasto(ia)) then

            do j=1,nu
              vec(j) = Hdense(i,j)
            enddo
            call cholsl(Sdense,nu,maxnh,paux,vec,vec)

C  vec contains the elements of H' in reduced space.
C  now H' must be formed in sparse format -

c  global indes of orbital i of reduced space
            ii = listh(i,mu)
            do j=1,nu
c  global indes of orbital j of reduced space
              jj = listh(j,mu)
              do k = 1, numh(ii)
                if (listh(k,ii) .eq. jj) hbar(k,ii)=vec(j)
              enddo
            enddo

          endif
C -
        enddo
         
C ...

      enddo
C ............

C *** Compute smallest and largest eigenvalues (Lanczos Method) ***

      call lanc1(2,hbar,numh,listh,nhmax,nbasis,emin)
      call lanc1(1,hbar,numh,listh,nhmax,nbasis,emax)

      emean=0.5*(emax+emin)
      deltae=0.55*(emax-emin)

C *** Calculate Chemical Potential using the Projection Method of
C              Goedecker (PRB 51,9455 (95)). ***
C 
C rr(in,ip) stores the in-th element of the vector resulting from 
C the application of the ip-th Chebyshev polynomial to the in-th atomic
C orbital. This is all what is needed to calculate the number of 
C electrons.

C  First scale and shift the hamiltonian ...

      do j=1,nbasis
        do i=1,numh(j)
          if (listh(i,j).eq.j) then
            hbar(i,j)=(hbar(i,j)-emean)/deltae
          else
            hbar(i,j)=hbar(i,j)/deltae
          endif
        enddo
      enddo     

C ...

C Calculate maximum length in unit cell ...
      rmax = 0.0d0
      do i = -1,1
        do j = -1,1
          do k = -1,1
            ri(1) = i*cell(1,1) + j*cell(2,1) + k*cell(3,1)
            ri(2) = i*cell(1,2) + j*cell(2,2) + k*cell(3,2)
            ri(3) = i*cell(1,3) + j*cell(2,3) + k*cell(3,3)
            rrmod = sqrt( ri(1)**2 + ri(2)**2 + ri(3)**2 )
            if (rrmod .gt. rmax) rmax = rrmod
          enddo
        enddo
      enddo
C ...

C initialize routine for neighbour search
      if (2.*rcoor .lt. rmax) then
        nna = nnmax
        call neighb(cell,rcoor,natoms,xa,0,0,nna,jan,xij,r2ij)
      endif

C initialize control vectors to zero 
      do i = 1,nbasis
        listvt(i) = 0
      enddo

C Loop over atoms ...............
      do ia = 1,natoms

        if (2.*rcoor .lt. rmax) then
C  look for neighbors of atom ia
          nna = nnmax
          call neighb(cell,rcoor,natoms,xa,ia,0,nna,jan,xij,r2ij)
          if (nna .gt. nnmax) then
            write(6,*) 'chempot: ERROR: nnmax = ',nnmax
            write(6,*) '                should be al least',nna
            write(6,*) '                Recompile'
            stop
          endif
        else
          nna = natoms
          do jj = 1,natoms
            jan(jj) = jj
          enddo
        endif

C build structure of sparse vector v ...

c clear list ot atoms considered within loc. range ...
        do jj=1,nna
          indexloc(jj) = 0
        enddo
        numloc = 0
c ...


        numv=0
        do 30 j = 1,nna
          ja = jan(j)

c check if ja has already been included in current vector ...
          do jj = 1,numloc
            if (ja .eq. indexloc(jj)) goto 30
          enddo
          numloc = numloc+1
          indexloc(numloc)=ja
c ...


          do jorb = 1,lasto(ja) - lasto(ja-1)
            nu = jorb + lasto(ja-1)
            numv=numv+1
            if (numv .gt. nvmax) then
              write(6,*) 'chempot: ERROR: nvmax = ',nvmax
              write(6,*) '                should be al least',numv
              write(6,*) '                Recompile'
              stop
            endif
            listv(numv) = nu
            listvt(nu) = numv
          enddo
30      continue
C ...

c       number of orbitals of atom ia
        norb = lasto(ia) - lasto(ia-1)

c     loop over orbitals of atom ia ...
        do iorb = 1,norb
          mu = iorb + lasto(ia-1)

          do 35 m=1,numv
            v(m,1)=0.
            v(m,2)=0.
 35       continue  

          imu = listvt(mu)

          v(imu,1)=1.0
          rr(mu,0)=v(imu,1)

          do j=1,numv
            ji = listv(j)
            numhp(ji)=0
            do i=1,numh(ji)
              m=listh(i,ji)
              jj = listvt(m)
              if (jj .ne. 0) then
                numhp(ji)=numhp(ji)+1
                listhp(numhp(ji),ji)=jj
                listhpp(numhp(ji),ji)=i
                v(j,2)=v(j,2)+hbar(i,ji)*v(jj,1)
              endif
            enddo
          enddo

          rr(mu,1)=v(imu,2)

          do 40 k=2,p-1

	    do j=1,numv
	      v(j,3)=-v(j,1)
              ji = listv(j)
              do i=1,numhp(ji)
                jj=listhp(i,ji)
                ii=listhpp(i,ji)
                v(j,3)=v(j,3)+2.*hbar(ii,ji)*v(jj,2)
              enddo
            enddo


            rr(mu,k)=v(imu,3)

	    j=0
            do i=1,numv
              v(i,1)=v(i,2)
                v(i,2)=v(i,3)
            enddo

 40       continue  


        enddo
c ...

c reset control vectors to cero
        do i=1,numv
          imu = listv(i)
          listvt(imu) = 0
        enddo

      enddo
C ............

      tol=0.0001
C Calculate chemical potential as the root of Nel - Tr(rho) = 0 ...

C Inverse temperature (in units of energy scaled so that the spectrum
C  lays between (-1,+1)
      betap = beta * deltae
C	write(6,*) 'betap',betap
C	write(6,*) 'rcoor',rcoor
C ............

C      do p=1,pmax
       chpotsh=zbrent(numb,mu1,mu2,tol)
         
C ...

C Shift and scale the result to absolute energy scale ...
       chpot=chpotsh*deltae+emean
C ...

C      write(20,*) p,chpot*13.6
C      enddo

C  CALCULATION OF THE GAP DISABLED; IT DOES NOT WORK PROPERLY

      return

CC Scale and shift the hamiltonian to absolute energy scale ...
C
C      do j=1,nbasis
C        do i=1,numh(j)
C          if (listh(i,j).eq.j) then
C            hbar(i,j)=deltae*hbar(i,j)+emean
C          else
C            hbar(i,j)=deltae*hbar(i,j)
C          endif
C        enddo
C      enddo
C
CC ...
C
C
CC  *** The gap is obtained by applying the Lanczos Method to
CC        the Folded Spectrum Method. See:
CC        Capaz-Koiler, J. Appl. Phys, 74, 5531 (93)
CC        Grosso et al, Nuovo Cimento D 15, 269 (93)
CC        Wang-Zunger, J. Chem. Phys. 100, 2394 (94) ***
C
C      eref1=chpot
C
CC Shift Hamiltonian...
C      do j=1,nbasis
C        do i=1,numh(j)
C          if (listh(i,j).eq.j) then
C            hbar(i,j)=hbar(i,j)-eref1
C          endif
C        enddo
C      enddo
CC ...
C
CC Solve (H-eref1)**2 by Lanczos...
C      call lanc2(hbar,numh,listh,nhmax,nbasis,eig1)
CC ...
C
C      eig1 = eig1 + eref1
C
C      delta=eig1-eref1
C      eref2=eref1-delta
C
C50    continue
C
CC Shift Hamiltonian ...
C      do j=1,nbasis
C        do i=1,numh(j)
C          if (listh(i,j).eq.j) then
C            hbar(i,j)=hbar(i,j)+eref1-eref2
C          endif
C        enddo
C      enddo
CC ...
C
CC Solve (H-eref2)**2 by Lanczos ...
C      call lanc2(hbar,numh,listh,nhmax,nbasis,eig2)
CC ...
C
C      eig2=eig2+eref2
C      gap=abs(eig1-eig2)
C
CC Check that levels are above and below the Fermi Level ...
C      if ((eig1 .gt. chpot .and. eig2 .gt. chpot) .or.
C     .    (eig1 .lt. chpot .and. eig2 .lt. chpot)) then
C        eref1=eref2
C        eref2=eref2-delta
C        goto 50
C      endif
CC ...
C
CC Convert to absolute energy scale ...
Cc      eig1 = eig1*deltae + emean
Cc      eig2 = eig2*deltae + emean
Cc      gap = gap*deltae
CC ...
C
CC Assign HOMO and LUMO ...
C      if (eig1. gt. eig2) then
C        homo=eig2
C        lumo=eig1
C      else
C        homo=eig1
C        lumo=eig2
C      endif
CC ...
C
C      return
      end        



      function numb(mu)
C **********************************************************************
C This function calculates the difference between the true number of
C electrons of the system, and the output number of electrons for a
C given Chemical Potential mu.
C
C Written by Maider Machado and P.Ordejon, June'98
C **********************************************************************
      implicit none
      include 'ordern.h'
      real*8 mu,Ne,numb,qtot
      integer ipmax,k,n,nbasis,p
      parameter(ipmax=200)
      real*8 rr(maxo,0:ipmax),c(ipmax)
      common /vect/rr,c,qtot,p,nbasis

      call chebfd(p,mu,c) 
c      write(16,*) qtot,nbasis
c      do ix=1,1001
c      x = -1. + 2.*(ix-1)/1000.
      Ne=0.

      do n=1,nbasis
        Ne=Ne+0.5*c(1)*rr(n,0)+c(2)*rr(n,1)
        do k=1,p-2
           Ne=Ne+c(k+2)*rr(n,k+1) 
        enddo
      enddo

      numb=qtot-2.*Ne

c      Ne=Ne+0.5*c(1)+c(2)*x
c      txm1 = 1
c      tx   = x
c      do k=1,p-2
c         txp1 = 2*x*tx - txm1
c         Ne=Ne+c(k+2)*txp1
c         txm1 = tx
c         tx   = txp1
c      enddo
c      write(7,*) x,Ne
c      enddo


      return
      end
     

      subroutine chebfd(n,mu,c)
C ***********************************************************************
C Calculates the coefficients of the Chebyshev polynomials
C expansion of the Fermi-Dirac function.
C Ref: W.H.Press et al. Numerical Recipes, Cambridge Univ. Press.
C
C Adapted by Maider Machado and P. Ordejon,  June'98
C ****************************** INPUT **********************************
C integer n                    : order of the expansion
C real*8 mu                    : chemical potential for F-D function
C ****************************** OUTPUT *********************************
C real*8 c(n)                  : expansion coefficients
C ***********************************************************************
  
      INTEGER N    
      INTEGER NMAX      
      REAL*8 B,A
      REAL*8 BMA,BPA
      REAL*8 C(N),PI, BETA, MU
      PARAMETER(NMAX=200,PI=3.141592653589793D0)
      PARAMETER (A=-1.,B=1.) 
      COMMON /BETA/BETA
      INTEGER J,K
      REAL*8 FAC, Y, F(NMAX)
      REAL*8 SUM
      BMA=0.5*(B-A)
      BPA=0.5*(B+A)
      DO K=1,N
        Y= COS(PI*(K-0.5)/N)
        IF (BETA*(Y*BMA+BPA-MU) .GT. 50.)  THEN
          F(K) = 0.0d0
        ELSE IF (BETA*(Y*BMA+BPA-MU) .LT. -50.) THEN
          F(K) = 1.0d0
        ELSE
          F(K)=1.0D0/(EXP(BETA*(Y*BMA+BPA-MU))+1.0D0)     
        ENDIF
      ENDDO
      FAC=2./N
      DO 30 J=1,N
        SUM=0.D0
        DO 20 K=1,N
	  SUM=SUM+F(K)*COS(PI*(J-1)*(K-0.5D0)/N)
  20    CONTINUE
        C(J) =FAC*SUM
  30  CONTINUE  
      RETURN
      END



      FUNCTION ZBRENT(FUNC,X1,X2,TOL)
c  ***************************************************
c Finds the root of a function FUNC known to lie
c between X1 and X2. The root, returned as zbrent,
c will be refined until its accuracy is tol.
c
C Converted to double precision from same routine in
C "Numerical Recipes", W.Press et al, Cambridge U.P.
c  ***************************************************     
      IMPLICIT NONE
      INTEGER ITMAX
      REAL*8 ZBRENT,TOL,X1,X2,FUNC,EPS
      EXTERNAL FUNC
      PARAMETER (ITMAX=300,EPS=3.E-4)
      INTEGER ITER
      REAL*8 A,B,C,D,E,FA,FB,FC,P,Q,R,
     *S,TOL1,XM

c      integer ix

c      do ix=1,1001
c        a=(x2-x1)*(ix-1)/1000. + x1
c        fa=func(a)
c	write(6,*) a,func(a)
c      enddo

      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.) STOP 'ZBRENT: Root must be bracketed'
      C=B
      D=B-A
      E=D
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
	FB=FUNC(B)
11    CONTINUE
      STOP 'ZBRENT exceeding maximum iterations.'
C      ZBRENT=B
C      RETURN
      END
 


      subroutine choldc(a,n,np,p)
C Cholesky decompositin of symmetric matrix
C Converted to double precision from same routine in
C "Numerical Recipes", W.Press et al, Cambridge U.P.
      implicit none
      integer n,np
      real*8 a(np,np),p(n)
      integer i,j,k
      real*8 sum

      do i=1,n
        do j=i,n
          sum=a(i,j)
          do k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
          enddo
          if (i .eq. j) then
            if (sum .le. 0.) stop 'choldc failed'
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
        enddo
      enddo
      return
      end

      subroutine cholsl(a,n,np,p,b,x)
C Solves linear system for a symmetric matrix 
C (in Cholesky form, as output from output of choldc)
C Converted to double precision from same routine in
C "Numerical Recipes", W.Press et al, Cambridge U.P.
      implicit none
      integer n,np
      real*8 a(np,np),b(n),p(n),x(n)
      integer i,k
      real*8 sum

      do i=1,n
        sum=b(i)
        do k=i-1,1,-1
          sum=sum-a(i,k)*x(k)
        enddo
        x(i)=sum/p(i)
      enddo
      do i=n,1,-1
        sum=x(i)
        do k=i+1,n
          sum=sum-a(k,i)*x(k)
        enddo
        x(i)=sum/p(i)
      enddo
      return
      end




      subroutine lanc1(opt,hbar,numh,listh,
     .                   nhmax,nbasis,ener)
C *********************************************************************
C Routine to calculate the mimimum or maximum eigenvalues of
C a given sparse Hamiltonian, by (2nd order) the Lanczos Method.
C
C Written by Maider Machado and P.Ordejon,  June'98
C ******************************** INPUT ******************************
C integer opt                     : 1 = compute minimun eigenval.
C                                   2 = compute maximum eigenval.
C real*8 hbar(maxnh,maxo)         : Sparse Hamiltonian
C integer numh(nhmax)             : control vector of hbar
C integer listh(nhmax,nbasis)     : control vector of hbar
C integer nhmax                   : Maximum number of nonzero elements of
C                                   each column of the hamiltonian
C integer nbasis                  : number of basis orbitals
C ******************************* OUTPUT *******************************
C real*8 ener                     : Eigenvalue
C **********************************************************************

      implicit none

      include 'ordern.h'

      integer
     .  nhmax, nbasis, opt

      integer 
     .  listh(nhmax,nbasis), numh(nbasis)

      real*8
     .  ener, hbar(maxnh,maxo)

C  Internal variables ...

      integer
     .  itmax
      parameter (itmax=200)

      real
     .  tol
      parameter (tol=0.0001d0)

      integer
     .  i, ii, j, k

      real*8
     .  a0, b1, b2, c1, diff, ecomp, eivec(2), 
     .  hv0(maxo), mod, norm, ran3, u0(maxo), v0(maxo), wr(2)

      external
     .   ran3

C ...

C  An unlikely number ...
      ecomp=4321.0987d0
C ...

C  Generate an initial normalized random vector ........
      mod=0.
      do i=1,nbasis
        u0(i)=2.*ran3(-i)-1.
        mod=mod+u0(i)**2 
      enddo
      mod=sqrt(mod)
      do i=1,nbasis
        u0(i)=u0(i)/mod
      enddo
C ...

      do i=1,nbasis
        v0(i)=0.
      enddo

C  Lanczos loop ................
      do k=1,itmax
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(i,j)
              v0(ii)=v0(ii)+hbar(i,j)*u0(j)
           enddo
        enddo
        a0=0.
        do i=1,nbasis
           a0=a0+u0(i)*v0(i)
        enddo
        norm=0.
        do i=1,nbasis
           v0(i)=v0(i)-a0*u0(i)
           norm=norm+v0(i)**2
        enddo
        b1=sqrt(norm)

        do i=1,nbasis
           hv0(i)=0.0
        enddo

        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(i,j)
              hv0(ii)=hv0(ii)+hbar(i,j)*v0(j)
           enddo
        enddo
        b2=0.
        do j=1,nbasis
           b2=b2+u0(j)*hv0(j)/b1
        enddo
        c1=0.
        do i=1,nbasis
           c1=c1+v0(i)*hv0(i)/norm
        enddo


C eigenvalues and eigenvectors ...

        wr(1) = 0.5d0*((a0+c1) 
     .          + sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
        wr(2) = 0.5d0*((a0+c1) 
     .          - sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))

        eivec(1)=1/sqrt(1+((a0+b1-wr(opt))/(b2+c1-wr(opt)))**2)
        eivec(2)=-eivec(1)*(a0+b1-wr(opt))/(b2+c1-wr(opt))

        ener=wr(opt)
C ...

 
        norm=0.
        diff=abs((wr(opt)-ecomp)/wr(opt))
        if (diff.gt.tol) then
          do i=1,nbasis 
             u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
             norm=norm+u0(i)**2
             v0(i)=0.
          enddo
          do j=1,nbasis
            u0(j)=u0(j)/sqrt(norm)
          enddo
          ecomp=wr(opt)
        else 
          goto 20
        endif
      enddo
C .................


CC  Lanczos loop ................
C      do k=1,itmax
C        do j=1,nbasis
C           do i=1,numh(j)
C              ii=listh(i,j)
C              v0(ii)=v0(ii)+hbar(i,j)*u0(j)
C           enddo
C        enddo
C        a0=0.
C        do i=1,nbasis
C           a0=a0+u0(i)*v0(i)
C        enddo
C        norm=0.
C        do i=1,nbasis
C           v0(i)=v0(i)-a0*u0(i)
C           norm=norm+v0(i)**2
C        enddo
C        b1=sqrt(norm)
C
C        do i=1,nbasis
C           hv0(i)=0.0
C        enddo
C
C        do j=1,nbasis
C           do i=1,numh(j)
C              ii=listh(i,j)
C              hv0(ii)=hv0(ii)+hbar(i,j)*v0(j)
C           enddo
C        enddo
C        b2=0.
C        do j=1,nbasis
C           b2=b2+u0(j)*hv0(j)/b1
C        enddo
C        c1=0.
C        do i=1,nbasis
C           c1=c1+v0(i)*hv0(i)/norm
C        enddo
C
C
CC eigenvalues and eigenvectors ...
C
C        wr(1) = 0.5d0*((a0+c1) 
C     .          + sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
C        wr(2) = 0.5d0*((a0+c1) 
C     .          - sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
C
C        eivec(1)=1/sqrt(1+((a0+b1-wr(opt))/(b2+c1-wr(opt)))**2)
C        eivec(2)=-eivec(1)*(a0+b1-wr(opt))/(b2+c1-wr(opt))
C
C        ener=wr(opt)
CC ...
C
C 
C        norm=0.
C        diff=abs((wr(opt)-ecomp)/wr(opt))
C        if (diff.gt.tol) then
C          do i=1,nbasis 
C             u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
C             norm=norm+u0(i)**2
C             v0(i)=0.
C          enddo
C          do j=1,nbasis
C            u0(j)=u0(j)/sqrt(norm)
C          enddo
C          ecomp=wr(opt)
C        else 
C          goto 20
C        endif
C      enddo
CC .................

      write(6,*) 'WARNING: lanc1 not converged after ',itmax,
     .           ' iterations'

20    return
      end





      subroutine lanc2(hbar,numh,listh,nhmax,nbasis,ener)
C *********************************************************************
C Routine to calculate the eigenvalue closest to cero for a given
C sparse Hamiltonian H, using the the Folded Spectrum Method
C (by (2nd order) the Lanczos Method).
C
C Written by Maider Machado and P.Ordejon,  June'98
C ******************************** INPUT ******************************
C real*8 hbar(maxnh,maxo)         : Sparse Hamiltonian
C integer numh(nhmax)             : control vector of hbar
C integer listh(nhmax,nbasis)     : control vector of hbar
C integer nhmax                   : Maximum number of nonzero elements of
C                                   each column of the hamiltonian
C integer nbasis                  : number of basis orbitals
C ******************************* OUTPUT *******************************
C real*8 ener                     : Eigenvalue of H
C **********************************************************************

      implicit none

      include 'ordern.h'

      integer
     .  nhmax,nbasis

      integer 
     .  listh(nhmax,nbasis), numh(nbasis)

      real*8
     .  ener, hbar(maxnh,maxo)

C  Internal variables ...

      integer
     .  itmax
      parameter (itmax=500)

      real
     .  tol
      parameter (tol=0.000001d0)

      integer
     .  i, ii, ij, j, k

      real*8
     .  a0, b1, b2, c1, diff, ecomp, eivec(2), 
     .  hv0(maxo), mod, norm, ran3, u0(maxo), v0(maxo), 
     .  v00(maxo), wr

      external
     .   ran3

C ...


C  An unlikely number ...
      ecomp=4321.0987d0
C ...

C  Generate an initial normalized random vector ........
      mod=0.
      do i=1,nbasis
        u0(i)=2.*ran3(-i)-1.
        mod=mod+u0(i)**2 
      enddo
      mod=sqrt(mod)
      do i=1,nbasis
        u0(i)=u0(i)/mod
      enddo
C ...

      do i=1,nbasis
        v0(i)=0.
        v00(i)=0.
      enddo

C  Lanczos loop ................
      do k=1,itmax
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(i,j)
              v00(j)=v00(j)+hbar(i,j)*u0(ii)
           enddo
        enddo
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(i,j)
              v0(j)=v0(j)+hbar(i,j)*v00(ii)
           enddo
        enddo
        a0=0.
        do i=1,nbasis
           a0=a0+u0(i)*v0(i)
        enddo
        norm=0.
        do i=1,nbasis
           v0(i)=v0(i)-a0*u0(i)
           norm=norm+v0(i)**2
        enddo
        b1=sqrt(norm)

        do i=1,nbasis
           hv0(i)=0.0
	   v00(i)=0.0
        enddo

        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(i,j)
              v00(j)=v00(j)+hbar(i,j)*v0(ii)
           enddo
        enddo
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(i,j)
              hv0(j)=hv0(j)+hbar(i,j)*v00(ii)
           enddo
        enddo
        b2=0.
        do j=1,nbasis
           b2=b2+u0(j)*hv0(j)/b1
        enddo
        c1=0.
        do i=1,nbasis
           c1=c1+v0(i)*hv0(i)/norm
        enddo


c  minimum eigenvalue ...
        wr = 0.5d0*((a0+c1) 
     .                 - sqrt((a0+c1)**2 - 4.0d0*(a0*c1-b1*b2)))
c ...

c eigenvector ...
        eivec(1)=1/sqrt(1+((a0+b1-wr)/(b2+c1-wr))**2)
        eivec(2)=-eivec(1)*(a0+b1-wr)/(b2+c1-wr)
C ...

 
        norm=0.

        do i=1,nbasis 
          u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
          norm=norm+u0(i)**2
          v0(i)=0.
          v00(i)=0.
        enddo
        do j=1,nbasis
          u0(j)=u0(j)/sqrt(norm)
        enddo

        diff=abs(wr-ecomp)
        ecomp=wr
        if (diff.lt.tol) goto 20
      enddo
C .................

      write(6,*) 'WARNING: lanc2 not converged after ',itmax,
     .           ' iterations'

20    continue

C Calculate eigehvalue of H ...
      ener=0.0
      do i=1,nbasis
        v0(i)=0.0
      enddo
      do i=1,nbasis
        do j=1,numh(i)
          ij=listh(j,i)
          v0(i)=v0(i)+hbar(j,i)*u0(ij)
        enddo
      enddo
      do i=1,nbasis
        ener=ener+u0(i)*v0(i)
      enddo
      
      return
      end




