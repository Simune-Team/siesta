      subroutine chempot(h,s,numh,listhptr,listh,rcoor,pmax,beta,
     .                   lasto,cell,xa,enum,nbasis,nbasisloc,
     .                   natoms,maxnh,chpot,emax,emin)

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
C       and gap is disabled because it doesnot work properly.
C
C NOTE: Parallelism has only been introduced to the extent of
C       handling the distribution of H and S. Needs more
C       extensive work to distribute local data to save
C       memory for large problems.
C
C Written by Maider Machado and P.Ordejon, June'98 
C ***************************** INPUT ***************************** 
C real*8 h(maxnh)                 : Hamiltonian in sparse form
C real*8 s(maxnh)                 : Overlap in sparse form
C integer numh(nbasisloc)         : Control vector of sparse hamilt.
C integer listhptr(nbasisloc)     : Control vector of sparse hamilt.
C integer listh(maxnh)            : Control vector of sparse hamilt.
C real*8 rcoor                    : Cutoff range for the projection 
C                                   vectors
C integer pmax                    : Maximum order of Chebishev expansion
C real*8 beta                     : Inverse Temperature for Chebi expansion
C integer lasto(0:natoms)         : Index vector of last orbital of 
C                                   each atom
C real*8 cell(3,3)                : Lattice vectors
C real*8 xa(3,natoms)             : Atomic positions
C real*8 enum                     : Total number of electrons
C integer nbasis                  : Number of basis orbitals (global)
C integer nbasisloc               : Number of basis orbitals (local)
C integer natoms                  : Number of atoms
C integer maxnh                   : Maximum number of non-zero elements
C                                   of sparse Hamiltonian
C **************************** OUTPUT *******************************
C real*8 chpot                    : Chemical potential
C real*8 emax                     : Maximum eigenvalue
C real*8 emin                     : Minimum eigenvalue
C *************************** DISABLED ******************************
C real*8 gap                      : Energy gap
C real*8 homo                     : Highest occ. molecular orbital energy
C real*8 lumo                     : Lowest unocc. molecular orbital energy
C *********************************************************************

C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi
#endif
      use numbvect

      implicit none

      integer 
     .  natoms, nbasis, nbasisloc, maxnh, pmax

      integer
     .  lasto(0:natoms), listh(maxnh), numh(nbasisloc), 
     .  listhptr(nbasisloc)

      real*8
     .  beta, cell(3,3), chpot, emin, emax, enum, 
     .  h(maxnh), rcoor, s(maxnh), xa(3,natoms)
C    .  homo, lumo, gap

C Internal variables...

      integer
     . ind, ipmax, nhmax, Node, Nodes, nnmaxnew, nvmaxnew
#ifdef MPI
      integer
     . MPIerror
#endif

C nnmax = maximum number of neighbors atoms within rcoor
C nvmax = maximum number of non-zero elements within a sparse vector v
      integer, save ::
     .  nnmax, nvmax

C ipmax = maximum order the Chebyshev polynomial exp. of the Fermi operator
      parameter(ipmax = 200)

      real*8 
     .  mu1,mu2
      parameter(mu1=-1.2d0,mu2=1.2d0)

      integer
     .  i, ia, ii, iil, imu, iorb, j, ja, ji, jj, jorb, 
     .  k, m, mu, nna, norb, nu, num, numloc, numv

      integer, dimension(:), allocatable, save ::
     .  listvt, numhp, indexloc, jan, ibuffer, listv

      real*8
     .  betap, chpotsh, deltae, emean,
     .  numb, rmax, rrmod, ri(3), tol,
     .  zbrent
C    .  delta, eref2, eref1, eig2, eig1

      integer, dimension(:,:), allocatable, save ::
     .  listhp, listhpp

      real*8, dimension(:), allocatable, save ::
     .  paux, vec

      real*8, dimension(:,:), allocatable, save ::
     .  hbar, Hdense, Sdense

      real*8, dimension(:), allocatable, save ::
     . r2ij, dpbuffer1

      real*8, dimension(:,:), allocatable, save ::
     . xij, v, dpbuffer2

C ...
C common to pass beta to chebfd ...
      common /beta/betap
C ...
      external 
     .  numb, zbrent, memory

      data nnmax / 1000 /
      data nvmax / 1000 /
C ...

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C Assign information for Chebyshev expansion .................
C Total charge
      qtot = enum
      nb = nbasis
      p = pmax

C Find size of maxnh
      nhmax = 0
      do i = 1,nbasisloc
        nhmax = max(nhmax,numh(i))
      enddo

C Allocate arrays that depend on maxnh
      allocate(listhp(nhmax,nbasis))
      call memory('A','I',nhmax*nbasis,'chempot')
      allocate(listhpp(nhmax,nbasis))
      call memory('A','I',nhmax*nbasis,'chempot')
      allocate(paux(nhmax))
      call memory('A','D',nhmax,'chempot')
      allocate(vec(nhmax))
      call memory('A','D',nhmax,'chempot')
      allocate(hbar(nhmax,nbasis))
      call memory('A','D',nhmax*nbasis,'chempot')
      allocate(Hdense(nhmax,nhmax))
      call memory('A','D',nhmax*nhmax,'chempot')
      allocate(Sdense(nhmax,nhmax))
      call memory('A','D',nhmax*nhmax,'chempot')

C Allocate arrays that depend on nbasis
      allocate(listvt(nbasis))
      call memory('A','I',nbasis,'chempot')
      allocate(numhp(nbasis))
      call memory('A','I',nbasis,'chempot')
      allocate(rr(nbasis,0:ipmax))
      call memory('A','D',nbasis*(ipmax+1),'chempot')

C Allocate arrays that depend on ipmax
      allocate(c(ipmax))
      call memory('A','D',ipmax,'chempot')

C Allocate arrays that depend on nnmax
      allocate(indexloc(nnmax))
      call memory('A','I',nnmax,'chempot')
      allocate(jan(nnmax))
      call memory('A','I',nnmax,'chempot')
      allocate(r2ij(nnmax))
      call memory('A','D',nnmax,'chempot')
      allocate(xij(3,nnmax))
      call memory('A','D',3*nnmax,'chempot')

C Allocate arrays that depend on nvmax
      allocate(listv(nvmax))
      call memory('A','I',nvmax,'chempot')
      allocate(v(nvmax,3))
      call memory('A','D',3*nvmax,'chempot')

C *** Calculate H' = S(-1)*H  (See Gibson et al, PRB 47, 9229 (92)) ***
C This is done by Choleski decomposition of S, and solving the
C linear system S H' = H, in the subspace of orbitals which overlap
C with those of a given atom. The advantage is that the Choleski
C decomposition can be done only once for each atom, and use
C the result for all the orbitals.
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
          call GlobalToLocalOrb(i,Node,Nodes,ii)
          if (numh(ii) .gt. num) then
            mu=ii
            num=numh(ii)
          endif
        enddo
        if (mu .eq. 0) then
          if (Node.eq.0) then
            write(6,*) 'chempot: ERROR: zero neighbors for orbital ',i
          endif
          stop
        endif
        nu=numh(mu) 
C -

C Initialize overlap and interaction in reduced space -
        do i=1,nu
          do j=1,nu
            Sdense(i,j)=0.0d0
            Hdense(i,j)=0.0d0
          enddo
        enddo
C -

C Construct S and H in reduced space  -
C  loop over orbitals ii in reduced space
        do i=1,nu
          ii = listh(listhptr(mu)+i)
          call GlobalToLocalOrb(ii,Node,Nodes,iil)
          if (iil.gt.0) then
C  loop over orbitals jj in reduced space
            do j=i,nu
              jj = listh(listhptr(mu)+j)
C  see if orbitals ii and jj interact
              do k = 1,numh(iil)
                ind = listhptr(iil)+k
                if (listh(ind) .eq. jj) then
                  Sdense(i,j)=s(ind)
                  Sdense(j,i)=Sdense(i,j)
                  Hdense(i,j)=h(ind)
                  Hdense(j,i)=Hdense(i,j)
                endif
              enddo
            enddo
          endif
        enddo

#ifdef MPI
C  Globalise Hdense and Sdense for now
        allocate(dpbuffer1(nu))
        call memory('A','D',nu,'chempot')
        do i=1,nu
          call MPI_AllReduce(Sdense(1,i),dpbuffer1,nu,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
          do j=1,nu
            Sdense(j,i) = dpbuffer1(j)
          enddo
          call MPI_AllReduce(Hdense(1,i),dpbuffer1,nu,DAT_double,
     .      MPI_sum,MPI_Comm_World,MPIerror)
          do j=1,nu
            Hdense(j,i) = dpbuffer1(j)
          enddo
        enddo
        call memory('D','D',size(dpbuffer1),'chempot')
        deallocate(dpbuffer1)
#endif

C  Cholesky factorization:
        call choldc(Sdense,nu,nhmax,paux)
C ...

C loop over orbitals of atom ia ....
C 
        do i=1,nu
          ii=listh(listhptr(mu)+i)
C  check if ii is in atom ia
          if (ii .ge. (lasto(ia-1)+1) .and. ii .le. lasto(ia)) then

            do j=1,nu
              vec(j) = Hdense(i,j)
            enddo
            call cholsl(Sdense,nu,nhmax,paux,vec,vec)

C  vec contains the elements of H' in reduced space.
C  now H' must be formed in sparse format -

c  global indes of orbital i of reduced space
            ii = listh(listhptr(mu)+i)
            call GlobalToLocalOrb(ii,Node,Nodes,iil)
            do j=1,nu
c  global indes of orbital j of reduced space
              jj = listh(listhptr(mu)+j)
              do k = 1, numh(iil)
                if (listh(listhptr(iil)+k) .eq. jj) hbar(k,ii)=vec(j)
              enddo
            enddo

          endif
C -
        enddo
C ...

      enddo
C ............

C *** Compute smallest and largest eigenvalues (Lanczos Method) ***

      call lanc1(2,hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,
     .  emin,Node)
      call lanc1(1,hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,
     .  emax,Node)

      emean=0.5d0*(emax+emin)
      deltae=0.55d0*(emax-emin)

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
          if (listh(listhptr(j)+i).eq.j) then
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
        if (nna .gt. nnmax) then
          nnmaxnew = nna + nint(0.1*nna)
C
          allocate(ibuffer(nnmax))
          call memory('A','I',nnmax,'chempot')
          ibuffer(1:nnmax) = indexloc(1:nnmax)
          call memory('D','I',size(indexloc),'chempot')
          deallocate(indexloc)
          allocate(indexloc(nnmaxnew))
          call memory('A','I',nnmaxnew,'chempot')
          indexloc(1:nnmax) = ibuffer(1:nnmax)
          ibuffer(1:nnmax) = jan(1:nnmax)
          call memory('D','I',size(jan),'chempot')
          deallocate(jan)
          allocate(jan(nnmaxnew))
          call memory('A','I',nnmaxnew,'chempot')
          jan(1:nnmax) = ibuffer(1:nnmax)
          call memory('D','I',size(ibuffer),'chempot')
          deallocate(ibuffer)
C
          allocate(dpbuffer1(nnmax))
          call memory('A','D',nnmax,'chempot')
          dpbuffer1(1:nnmax) = r2ij(1:nnmax)
          call memory('D','D',size(r2ij),'chempot')
          deallocate(r2ij)
          allocate(r2ij(nnmaxnew))
          call memory('A','D',nnmaxnew,'chempot')
          r2ij(1:nnmax) = dpbuffer1(1:nnmax)
          call memory('D','D',size(dpbuffer1),'chempot')
          deallocate(dpbuffer1)
C
          allocate(dpbuffer2(3,nnmax))
          call memory('A','D',3*nnmax,'chempot')
          dpbuffer2(1:3,1:nnmax) = xij(1:3,1:nnmax)
          call memory('D','D',size(xij),'chempot')
          deallocate(xij)
          allocate(xij(3,nnmaxnew))
          call memory('A','D',3*nnmaxnew,'chempot')
          xij(1:3,1:nnmax) = dpbuffer2(1:3,1:nnmax)
          call memory('D','D',size(dpbuffer2),'chempot')
          deallocate(dpbuffer2)
C
          nnmax = nnmaxnew
        endif
      endif

C initialize control vectors to zero 
      listvt(1:nbasis) = 0

C Loop over atoms ...............
      do ia = 1,natoms

        if (2.*rcoor .lt. rmax) then
C  look for neighbors of atom ia
          nna = nnmax
          call neighb(cell,rcoor,natoms,xa,ia,0,nna,jan,xij,r2ij)
          if (nna .gt. nnmax) then
            nnmaxnew = nna + nint(0.1*nna)
C
            allocate(ibuffer(nnmax))
            call memory('A','I',nnmax,'chempot')
            ibuffer(1:nnmax) = indexloc(1:nnmax)
            call memory('D','I',size(indexloc),'chempot')
            deallocate(indexloc)
            allocate(indexloc(nnmaxnew))
            call memory('A','I',nnmaxnew,'chempot')
            indexloc(1:nnmax) = ibuffer(1:nnmax)
            ibuffer(1:nnmax) = jan(1:nnmax)
            call memory('D','I',size(jan),'chempot')
            deallocate(jan)
            allocate(jan(nnmaxnew))
            call memory('A','I',nnmaxnew,'chempot')
            jan(1:nnmax) = ibuffer(1:nnmax)
            call memory('D','I',size(ibuffer),'chempot')
            deallocate(ibuffer)
C
            allocate(dpbuffer1(nnmax))
            call memory('A','D',nnmax,'chempot')
            dpbuffer1(1:nnmax) = r2ij(1:nnmax)
            call memory('D','D',size(r2ij),'chempot')
            deallocate(r2ij)
            allocate(r2ij(nnmaxnew))
            call memory('A','D',nnmaxnew,'chempot')
            r2ij(1:nnmax) = dpbuffer1(1:nnmax)
            call memory('D','D',size(dpbuffer1),'chempot')
            deallocate(dpbuffer1)
C
            allocate(dpbuffer2(3,nnmax))
            call memory('A','D',3*nnmax,'chempot')
            dpbuffer2(1:3,1:nnmax) = xij(1:3,1:nnmax)
            call memory('D','D',size(xij),'chempot')
            deallocate(xij)
            allocate(xij(3,nnmaxnew))
            call memory('A','D',3*nnmaxnew,'chempot')
            xij(1:3,1:nnmax) = dpbuffer2(1:3,1:nnmax)
            call memory('D','D',size(dpbuffer2),'chempot')
            deallocate(dpbuffer2)
C
            nnmax = nnmaxnew
          endif
        else
          nna = natoms
          do jj = 1,natoms
            jan(jj) = jj
          enddo
        endif

C build structure of sparse vector v ...

c clear list ot atoms considered within loc. range ...
        indexloc(1:nna) = 0
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
              nvmaxnew = numv + nint(0.1*numv)
C
              allocate(ibuffer(nvmax))
              call memory('A','I',nvmax,'chempot')
              ibuffer(1:nvmax) = listv(1:nvmax)
              call memory('D','I',size(listv),'chempot')
              deallocate(listv)
              allocate(listv(nvmaxnew))
              call memory('A','I',nvmaxnew,'chempot')
              listv(1:nvmax) = ibuffer(1:nvmax)
              call memory('D','I',size(ibuffer),'chempot')
              deallocate(ibuffer)
C
              allocate(dpbuffer2(nvmax,3))
              call memory('A','D',3*nvmax,'chempot')
              dpbuffer2(1:nvmax,1:3) = v(1:nvmax,1:3)
              call memory('D','D',size(v),'chempot')
              deallocate(v)
              allocate(v(nvmaxnew,3))
              call memory('A','D',3*nvmaxnew,'chempot')
              v(1:nvmax,1:3) = dpbuffer2(1:nvmax,1:3)
              call memory('D','D',size(dpbuffer2),'chempot')
              deallocate(dpbuffer2)
C
              nvmax = nvmaxnew
            endif
            listv(numv) = nu
            listvt(nu) = numv
          enddo
30      continue

c number of orbitals of atom ia
        norb = lasto(ia) - lasto(ia-1)

c loop over orbitals of atom ia ...
        do iorb = 1,norb
          mu = iorb + lasto(ia-1)

          v(1:numv,1:2)=0.0d0

          imu = listvt(mu)

          v(imu,1)=1.0d0
          rr(mu,0)=v(imu,1)

          do j=1,numv
            ji = listv(j)
            numhp(ji)=0
            do i=1,numh(ji)
              m=listh(listhptr(ji)+i)
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

      tol=0.0001d0
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

C      write(20,*) p,chpot*13.6d0
C      enddo

C Deallocate local arrays 
      call memory('D','I',size(listhp),'chempot')
      deallocate(listhp)
      call memory('D','I',size(listhpp),'chempot')
      deallocate(listhpp)
      call memory('D','D',size(paux),'chempot')
      deallocate(paux)
      call memory('D','D',size(vec),'chempot')
      deallocate(vec)
      call memory('D','D',size(hbar),'chempot')
      deallocate(hbar)
      call memory('D','D',size(Hdense),'chempot')
      deallocate(Hdense)
      call memory('D','D',size(Sdense),'chempot')
      deallocate(Sdense)
      call memory('D','I',size(listvt),'chempot')
      deallocate(listvt)
      call memory('D','I',size(numhp),'chempot')
      deallocate(numhp)
      call memory('D','D',size(rr),'chempot')
      deallocate(rr)
      call memory('D','D',size(c),'chempot')
      deallocate(c)
      call memory('D','I',size(indexloc),'chempot')
      deallocate(indexloc)
      call memory('D','I',size(jan),'chempot')
      deallocate(jan)
      call memory('D','D',size(r2ij),'chempot')
      deallocate(r2ij)
      call memory('D','D',size(xij),'chempot')
      deallocate(xij)
      call memory('D','I',size(listv),'chempot')
      deallocate(listv)
      call memory('D','D',size(v),'chempot')
      deallocate(v)

C  CALCULATION OF THE GAP DISABLED; IT DOES NOT WORK PROPERLY

      return

CC Scale and shift the hamiltonian to absolute energy scale ...
C
C      do j=1,nbasis
C        do i=1,numh(j)
C          if (listh(listhptr(j)+i).eq.j) then
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
C          if (listh(listhptr(j)+i).eq.j) then
C            hbar(i,j)=hbar(i,j)-eref1
C          endif
C        enddo
C      enddo
CC ...
C
CC Solve (H-eref1)**2 by Lanczos...
C      call lanc2(hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,eig1,Node)
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
C          if (listh(listhptr(j)+i).eq.j) then
C            hbar(i,j)=hbar(i,j)+eref1-eref2
C          endif
C        enddo
C      enddo
CC ...
C
CC Solve (H-eref2)**2 by Lanczos ...
C      call lanc2(hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,eig2,Node)
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
      use numbvect

      implicit none

      real*8 mu,Ne,numb
      integer k,n

      call chebfd(p,mu,c) 
c      write(16,*) qtot,nb
c      do ix=1,1001
c      x = -1. + 2.*(ix-1)/1000.
      Ne=0.0d0

      do n=1,nb
        Ne=Ne+0.5*c(1)*rr(n,0)+c(2)*rr(n,1)
        do k=1,p-2
           Ne=Ne+c(k+2)*rr(n,k+1) 
        enddo
      enddo

      numb=qtot-2.0d0*Ne

c      Ne=Ne+0.5d0*c(1)+c(2)*x
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

      subroutine lanc1(opt,hbar,nhmax,numh,listhptr,listh,
     .                   maxnh,nbasis,ener,Node)
C *********************************************************************
C Routine to calculate the mimimum or maximum eigenvalues of
C a given sparse Hamiltonian, by (2nd order) the Lanczos Method.
C
C Written by Maider Machado and P.Ordejon,  June'98
C ******************************** INPUT ******************************
C integer opt                     : 1 = compute minimun eigenval.
C                                   2 = compute maximum eigenval.
C real*8 hbar(maxnh,nbasis)       : Sparse Hamiltonian
C integer nhmax                   : Lower dimension of hbar
C integer numh(nhmax)             : control vector of hbar
C integer listhptr(nhmax)         : control vector of hbar
C integer listh(maxnh)            : control vector of hbar
C integer maxnh                 : Maximum number of nonzero elements of
C                                   the hamiltonian
C integer nbasis                  : number of basis orbitals
C integer Node                    : local node number
C ******************************* OUTPUT *******************************
C real*8 ener                     : Eigenvalue
C **********************************************************************

      implicit none

      integer
     .  nhmax, nbasis, opt, Node, maxnh

      integer 
     .  listh(maxnh), numh(nbasis), listhptr(nbasis)

      real*8
     .  ener, hbar(nhmax,nbasis)

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
     .  mod, norm, ran3, wr(2)

      real*8, dimension(:), allocatable, save ::
     .  hv0, u0, v0

      external
     .   ran3

C ...

C  Allocate local memory
      allocate(hv0(nbasis))
      call memory('A','D',nbasis,'chempot')
      allocate(v0(nbasis))
      call memory('A','D',nbasis,'chempot')
      allocate(u0(nbasis))
      call memory('A','D',nbasis,'chempot')

C  An unlikely number ...
      ecomp=4321.0987d0
C ...

C  Generate an initial normalized random vector ........
      mod=0.0d0
      do i=1,nbasis
        u0(i)=2.0d0*ran3(-i)-1.
        mod=mod+u0(i)**2 
      enddo
      mod=sqrt(mod)
      do i=1,nbasis
        u0(i)=u0(i)/mod
      enddo
C ...

      v0(1:nbasis)=0.0d0

C  Lanczos loop ................
      do k=1,itmax
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v0(ii)=v0(ii)+hbar(i,j)*u0(j)
           enddo
        enddo
        a0=0.0d0
        do i=1,nbasis
           a0=a0+u0(i)*v0(i)
        enddo
        norm=0.0d0
        do i=1,nbasis
           v0(i)=v0(i)-a0*u0(i)
           norm=norm+v0(i)**2
        enddo
        b1=sqrt(norm)

        hv0(1:nbasis)=0.0d0

        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              hv0(ii)=hv0(ii)+hbar(i,j)*v0(j)
           enddo
        enddo
        b2=0.0d0
        do j=1,nbasis
           b2=b2+u0(j)*hv0(j)/b1
        enddo
        c1=0.0d0
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

        norm=0.0d0
        diff=abs((wr(opt)-ecomp)/wr(opt))
        if (diff.gt.tol) then
          do i=1,nbasis 
             u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
             norm=norm+u0(i)**2
             v0(i)=0.0d0
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
C              ii=listh(listhptr(j)+i)
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
C              ii=listh(listhptr(j)+i)
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

      if (Node.eq.0) then
        write(6,*) 'WARNING: lanc1 not converged after ',itmax,
     .           ' iterations'
      endif

20    continue

C  Deallocate local memory
      call memory('D','D',size(hv0),'chempot')
      deallocate(hv0)
      call memory('D','D',size(v0),'chempot')
      deallocate(v0)
      call memory('D','D',size(u0),'chempot')
      deallocate(u0)

      return
      end


      subroutine lanc2(hbar,nhmax,numh,listhptr,listh,maxnh,nbasis,
     .                 ener,Node)
C *********************************************************************
C Routine to calculate the eigenvalue closest to cero for a given
C sparse Hamiltonian H, using the the Folded Spectrum Method
C (by (2nd order) the Lanczos Method).
C
C Written by Maider Machado and P.Ordejon,  June'98
C ******************************** INPUT ******************************
C real*8 hbar(maxnh,nbasis)         : Sparse Hamiltonian
C integer nhmax                   : Lower dimension of hbar
C integer numh(nhmax)             : control vector of hbar
C integer listhptr(nhmax)         : control vector of hbar
C integer listh(maxnh)            : control vector of hbar
C integer maxnh                   : Maximum number of nonzero elements of
C                                   the hamiltonian
C integer nbasis                  : number of basis orbitals
C integer Node                    : local node number
C ******************************* OUTPUT *******************************
C real*8 ener                     : Eigenvalue of H
C **********************************************************************

      implicit none

      integer
     .  nhmax, nbasis, Node, maxnh

      integer 
     .  listh(maxnh), numh(nbasis), listhptr(nbasis)

      real*8
     .  ener, hbar(nhmax,nbasis)

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
     .  mod, norm, ran3, wr

      real*8, dimension(:), allocatable, save ::
     .  hv0, u0, v0, v00

      external
     .   ran3

C ...

C  Allocate local memory
      allocate(hv0(nbasis))
      call memory('A','D',nbasis,'chempot')
      allocate(u0(nbasis))
      call memory('A','D',nbasis,'chempot')
      allocate(v0(nbasis))
      call memory('A','D',nbasis,'chempot')
      allocate(v00(nbasis))
      call memory('A','D',nbasis,'chempot')

C  An unlikely number ...
      ecomp=4321.0987d0
C ...

C  Generate an initial normalized random vector ........
      mod=0.0d0
      do i=1,nbasis
        u0(i)=2.0d0*ran3(-i)-1.
        mod=mod+u0(i)**2 
      enddo
      mod=sqrt(mod)
      do i=1,nbasis
        u0(i)=u0(i)/mod
      enddo
C ...

      v0(1:nbasis)=0.0d0
      v00(1:nbasis)=0.0d0

C  Lanczos loop ................
      do k=1,itmax
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v00(j)=v00(j)+hbar(i,j)*u0(ii)
           enddo
        enddo
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v0(j)=v0(j)+hbar(i,j)*v00(ii)
           enddo
        enddo
        a0=0.0d0
        do i=1,nbasis
           a0=a0+u0(i)*v0(i)
        enddo
        norm=0.0d0
        do i=1,nbasis
           v0(i)=v0(i)-a0*u0(i)
           norm=norm+v0(i)**2
        enddo
        b1=sqrt(norm)

        hv0(1:nbasis)=0.0d0
	v00(1:nbasis)=0.0d0

        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              v00(j)=v00(j)+hbar(i,j)*v0(ii)
           enddo
        enddo
        do j=1,nbasis
           do i=1,numh(j)
              ii=listh(listhptr(j)+i)
              hv0(j)=hv0(j)+hbar(i,j)*v00(ii)
           enddo
        enddo
        b2=0.0d0
        do j=1,nbasis
           b2=b2+u0(j)*hv0(j)/b1
        enddo
        c1=0.0d0
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

        norm=0.0d0

        do i=1,nbasis 
          u0(i)=eivec(1)*u0(i)+eivec(2)*v0(i)/b1
          norm=norm+u0(i)**2
          v0(i)=0.0d0
          v00(i)=0.0d0
        enddo
        do j=1,nbasis
          u0(j)=u0(j)/sqrt(norm)
        enddo

        diff=abs(wr-ecomp)
        ecomp=wr
        if (diff.lt.tol) goto 20
      enddo
C .................

      if (Node.eq.0) then
        write(6,*) 'WARNING: lanc2 not converged after ',itmax,
     .           ' iterations'
      endif

20    continue

C Calculate eigehvalue of H ...
      ener=0.0d0
      v0(1:nbasis)=0.0d0
      do i=1,nbasis
        do j=1,numh(i)
          ij=listh(listhptr(i)+j)
          v0(i)=v0(i)+hbar(j,i)*u0(ij)
        enddo
      enddo
      do i=1,nbasis
        ener=ener+u0(i)*v0(i)
      enddo
      
C  Deallocate local memory
      call memory('D','D',size(hv0),'chempot')
      deallocate(hv0)
      call memory('D','D',size(u0),'chempot')
      deallocate(u0)
      call memory('D','D',size(v0),'chempot')
      deallocate(v0)
      call memory('D','D',size(v00),'chempot')
      deallocate(v00)

      return
      end
