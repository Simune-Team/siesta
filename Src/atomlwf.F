      subroutine cspa(ioptlwf,iopt,natoms,nbasis,nbasisloc,lasto,
     .        isa,qa,enum,rcoor,rh,cell,xa,maxnh,numh,listhptr,
     .        listh,maxnc,c,numc,listc,ncmax,nctmax,nbands,overflow)
C ******************************************************************************
C This subroutine builds the Localized Wave Functions, centered
C on ATOMS (searching for atoms within a cutoff radius rcoor), and 
C assigns a RANDOM initial guess to start the CG minimization.
C
C Criterion to build LWF's: 
C 1) Method of Kim et al: use more localized orbitals than 
C    occupied orbitals.
C    We assign the minimum number of orbitals so that there
C    is place for more electrons than those in the system;
C    for instance:
C      H:        1 LWF
C      C,Si:     3 LWF's
C      N:        3 LWF's
C      O:        4 LWF's
C      ...
C    It can be used with charged cells, although it has not been
C    tested in detail.
C 2) Method of Ordejon et al: number of localized orbitals 
C    equal to number of occupied orbitals. 
C    Only available when each atom has an EVEN number of
C    electrons. It can NOT be used with charged cells.
C
C Writen by P.Ordejon, 1993. 
C Re-writen by P.Ordejon, November'96.
C Corrected by P.Ordejon, April'97,  May'97  
C lmax, lmaxs and nzls erased from the input by DSP, Aug 1998.
C ******************************* INPUT ***************************************
C integer ioptlwf           : Build LWF's according to:
C                             0 = Read blindly from disk
C                             1 = Functional of Kim et al.
C                             2 = Functional of Ordejon-Mauri
C integer iopt              : 0 = Find structure of sparse C matrix and
C                                  build initial guess
C                             1 = Just find structure of C matrix
C integer natoms            : Number of atoms
C integer nbasis            : Number of basis orbitals
C integer nbasisloc         : Number of local basis orbitals
C integer lasto(0:natoms)   : Index of last orbital of each atom
C integer isa(natoms)       : Species index of each atom
C real*8 qa(natoms)         : Neutral atom charge
C real*8 enum               : Total number of electrons
C real*8 rcoor              : Cutoff radius of Localized Wave Functions
C real*8 rh                 : Maximum cutoff radius of Hamiltonian matrix
C real*8 cell(3,3)          : Supercell vectors
C real*8 xa(3,natoms)       : Atomic coordinates
C integer maxnc             : First dimension of C matrix, and maximum
C                             number of nonzero elements of each row of C
C ****************************** OUTPUT **************************************
C real*8 c(ncmax,nbasis)    : Initial guess for sparse coefficients 
C                             of LWF's  (only if iopt = 0)
C integer numc(nbasis)      : Control vetor of C matrix
C integer listc(ncmax,nbasis): Control vetor of C matrix
C integer ncmax             : True value for ncmax, 
C                             If ncmax it is too small, then
C                             c, numc and listc are NOT initialized!!
C integer nctmax            : Maximum number of nonzero elements
C                             of eaxh column of C
C integer nbands            : Number of LWF's
C ************************** INPUT AND OUTPUT ********************************
C logical overflow          : True: Maximum dimensions too small
C ****************************************************************************

C
C  Modules
C
      use precision
      use parallel
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer 
     .  iopt,ioptlwf,maxnc,natoms,nbands,nbasis,ncmax,nctmax,
     .  maxnh,nbasisloc

      integer 
     .  isa(natoms),lasto(0:natoms),
     .  numc(nbasis),listc(maxnc,nbasis),numh(nbasisloc),
     .  listhptr(nbasisloc),listh(maxnh)

      double precision
     .  cell(3,3),enum,qa(natoms),rcoor,rh,xa(3,natoms),
     .  c(maxnc,nbasisloc)
 
      logical
     .  overflow

C  Internal variables .......................................................
C  nnmax = maximum number of neighbour atoms within rcoor
      integer, save ::
     . nnmax

      integer
     . i,ia,imu,
     . in,index,indexa,indexb,indexi,indexj,iorb,iseed,
     . j,ja,jj,k,mu,nelectr,nm,nna,norb,nqtot,nu,
     . numloc,Node,Nodes,mul,nct

      integer, dimension(:), allocatable, save ::
     . ilist, iptrgl, indexloc, jan

      logical
     .  lfound
#ifdef MPI
      integer
     .  MPIerror, ibuf1(2), ibuf2(2), imin, ntmp, nn, ntot
      integer, dimension(:), allocatable, save ::
     .  itmp, numcl
      double precision
     .  rtmp
      logical
     .  ltmp
#endif

      double precision
     . cg,fact,qtot,r,rmax,rr(3),rrmod,
     . snor,tiny,ran3

      double precision, dimension(:), allocatable, save ::
     . r2ij

      double precision, dimension(:,:), allocatable, save ::
     . xij

      data nnmax / 1000 /
C ...........................................................................

      call timer('cspa',1)

      tiny = 1.d-10

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

C Check that iopt is correct ................................................
      if (iopt .ne. 0 .and. iopt .ne. 1) then
        if (Node.eq.0) then
          write(6,*) 'cspa: Wrong iopt in cspa'
        endif
        stop
      endif
C ..................

C Allocate local memory
      allocate(ilist(nbasis))
      call memory('A','I',nbasis,'cspa')
      allocate(iptrgl(nbasis))
      call memory('A','I',nbasis,'cspa')
   10 allocate(indexloc(nnmax))
      call memory('A','I',nnmax,'cspa')
      allocate(jan(nnmax))
      call memory('A','I',nnmax,'cspa')
      allocate(r2ij(nnmax))
      call memory('A','D',nnmax,'cspa')
      allocate(xij(3,nnmax))
      call memory('A','D',3*nnmax,'cspa')

C ..................

C Initialize some stuff ......................................................
      do mu = 1,nbasis
        numc(mu) = 0
      enddo

      if (iopt .eq. 0) then
        do mu = 1,nbasisloc
          do i = 1,maxnc
            c(i,mu) = 0.0d0
          enddo
        enddo
      endif

C Need to reset random number generate here to ensure that all Nodes
C are going to access the same sequence of random numbers.
      iseed = -1
      cg = ran3(iseed)

      iseed = 17

      ncmax   = 0
      nctmax  = 0
C ........................

C Calculate maximum length in unit cell ......................................
C determine maximum cell length
      rmax = 0.0d0
      do i = -1,1
        do j = -1,1
          do k = -1,1
            rr(1) = i*cell(1,1) + j*cell(2,1) + k*cell(3,1)
            rr(2) = i*cell(1,2) + j*cell(2,2) + k*cell(3,2)
            rr(3) = i*cell(1,3) + j*cell(2,3) + k*cell(3,3)
            rrmod = sqrt( rr(1)**2 + rr(2)**2 + rr(3)**2 )
            if (rrmod .gt. rmax) rmax = rrmod
          enddo
        enddo
      enddo
C ........................

C Check that there is an even number of electrons in the system ...............
      nqtot = nint(enum + tiny)
      if (abs(nqtot-enum+tiny) .gt. 1d-3.and.Node.eq.0) then
        write(6,*) 'cspa: WARNING: total charge non integer:',enum
      endif
      if (abs(2.0*dble(nqtot/2)-dble(nqtot)).gt.1.0d-3) then
        if (Node.eq.0) then
        write(6,*) 'cspa: ERROR: Wrong total charge; odd charge:',nqtot
        write(6,*) '      ERROR: Charge must be EVEN to use Order-N ',
     .    'option'
        endif
        stop
      endif
C ..................

C Check if system is charged..............
      qtot = 0.0d0
      do ia = 1,natoms
        qtot = qtot + qa(ia)
      enddo
      if (abs(qtot-enum) .gt. 1.0d-4) then
        if (ioptlwf .eq. 2) then
          if (Node.eq.0) then
          write(6,*) 'cspa: ERROR: Charged systems can not be done'
          write(6,*) '      ERROR: with the Ordejon-Mauri functional'
          endif
          stop
        endif
      endif
C ..................

C Build control vectors of sparse LWF 
C loop over the localized wave funcions (centered at atoms)..................

C Initialize routine for neighbour search
      nna = nnmax
      call neighb(cell,rcoor,natoms,xa,0,0,nna,jan,xij,r2ij)
      if (nna .gt. nnmax) then
        call memory('D','I',size(indexloc),'cspa')
        deallocate(indexloc)
        call memory('D','I',size(jan),'cspa')
        deallocate(jan)
        call memory('D','D',size(r2ij),'cspa')
        deallocate(r2ij)
        call memory('D','D',size(xij),'cspa')
        deallocate(xij)
        nnmax = nna + nint(0.1*nna)
        goto 10
      endif

      index = 0
      do ia = 1,natoms

C  determine how many LWF's depending on the atomic species 
C  (or the number of electrons)
        nelectr = qa (ia) + tiny
        if (abs(nelectr - qa(ia) + tiny) .gt. 1e-3) then
          if (Node.eq.0) then
            write(6,*) 'cspa: Wrong atomic charge for atom ',ia
            write(6,*) '      qa = ',qa(ia),'  Must be an integer!!'
          endif
          stop
        endif
        if (ioptlwf .eq. 1) then
          indexi = ( ( nelectr + 2 ) / 2 )
        else if (ioptlwf .eq. 2) then
          if ( (nelectr/2)*2 .ne. nelectr) then
            if (Node.eq.0) then
              write(6,*) 'cspa: Wrong Order-N functional option in cspa'
              write(6,*) '      You can only use the functional of'
              write(6,*) '      Ordejon-Mauri for atoms with an even'
              write(6,*) '      number of electrons.'
            endif
            stop
          endif
          indexi = ( ( nelectr ) / 2 )
        else
          if (Node.eq.0) then
            write(6,*) 'cspa: Wrong functional option in cspa'
          endif
          stop
        endif

        if (mod((ia-1),Nodes).eq.Node) then
C Do this atom on the local node

          if (2.*rcoor .lt. rmax) then
C Look for neighbours of atom ia
            nna = nnmax
            call neighb(cell,rcoor,natoms,xa,ia,0,nna,jan,xij,r2ij)
            if (nna .gt. nnmax) then
              call memory('D','I',size(indexloc),'cspa')
              deallocate(indexloc)
              call memory('D','I',size(jan),'cspa')
              deallocate(jan)
              call memory('D','D',size(r2ij),'cspa')
              deallocate(r2ij)
              call memory('D','D',size(xij),'cspa')
              deallocate(xij)
              nnmax = nna + nint(0.1*nna)
              goto 10
            endif
          else
            nna = natoms
            do  jj = 1,natoms
              jan(jj) = jj
            enddo
          endif

C loop over LWF's centered on atom ia
	  do indexb = 1,indexi
	    index = index + 1
  
c clear list of atoms considered within loc. range
            do indexa = 1, nna
              indexloc(indexa) = 0
            enddo
            numloc=0

C initialize stuff...
            nct = 0       
 
C loop over the neighbors of ia within rcoor
            do  jj = 1,nna
              ja = jan(jj)
	      norb = lasto(ja)-lasto(ja-1)

C check if ja has already been included in current lwf
              do indexa = 1,numloc
                if (ja .eq. indexloc(indexa)) goto 30
              enddo
              numloc = numloc + 1
              indexloc(numloc) = ja

C loop over orbitals of ja
              do iorb = 1,norb
	        mu = iorb + lasto(ja-1)
                nct = nct + 1
	        nm = numc(mu)
	        numc(mu) = nm + 1
                if (numc(mu) .gt. maxnc) overflow = .true.
                if (.not.overflow) listc(nm+1,mu) = index
              enddo

30          enddo

            nctmax = max ( nctmax , nct )

          enddo

        else

          index = index + indexi

        endif

      enddo

#ifdef MPI
C Global reduction of nct
      call MPI_AllReduce(nctmax,nct,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)
      nctmax = nct

C Global reduction of overflow
      call MPI_AllReduce(overflow,ltmp,1,MPI_logical,MPI_lor,
     .  MPI_Comm_World,MPIerror)
      overflow = ltmp

      allocate(numcl(nbasis))
      if (.not.overflow) then
C Collect together numc and listc
        numcl(1:nbasis) = numc(1:nbasis)
        do nn = 1,Nodes
          if (Node.eq.nn-1) then
            ilist(1:nbasis) = numcl(1:nbasis)
          endif
          call MPI_Bcast(ilist,nbasis,MPI_integer,nn-1,
     .      MPI_Comm_World,MPIerror)
          ntot = 0
          do mu = 1,nbasis
            ntot = ntot + ilist(mu)
          enddo
          allocate(itmp(ntot))
          if (Node.eq.nn-1) then
            ntot = 0
            do mu = 1,nbasis
              itmp(ntot+1:ntot+ilist(mu)) = listc(1:ilist(mu),mu)
              ntot = ntot + ilist(mu)
            enddo
          endif
          call MPI_Bcast(itmp,ntot,MPI_integer,nn-1,
     .      MPI_Comm_World,MPIerror)
          if (Node.ne.nn-1) then
            ntot = 0
            do mu = 1,nbasis
              listc(numc(mu)+1:numc(mu)+ilist(mu),mu) = 
     .          itmp(ntot+1:ntot+ilist(mu))
              numc(mu) = numc(mu) + ilist(mu)
              ntot = ntot + ilist(mu)
            enddo
          endif
          deallocate(itmp)
        enddo

C Sort listc into order
        do mu = 1,nbasis
          do nn = 1,numc(mu)
            i = nn
            imin = listc(nn,mu)
            do in = nn+1,numc(mu)
              if (listc(in,mu).lt.imin) then
                i = in
                imin = listc(in,mu)
              endif
            enddo
            listc(i,mu) = listc(nn,mu)
            listc(nn,mu) = imin
          enddo
        enddo

      else
C Collect together numc and listc
        numcl(1:nbasis) = numc(1:nbasis)
        do nn = 1,Nodes
          if (Node.eq.nn-1) then
            ilist(1:nbasis) = numcl(1:nbasis)
          endif
          call MPI_Bcast(ilist,nbasis,MPI_integer,nn-1,
     .      MPI_Comm_World,MPIerror)
          if (Node.ne.nn-1) then
            do mu = 1,nbasis
              numc(mu) = numc(mu) + ilist(mu)
            enddo
          endif
        enddo
      endif
      deallocate(numcl)
#endif

      do mu = 1,nbasis
        ncmax  = max ( ncmax  , numc(mu) )
      enddo

      nbands = index

      if (iopt .eq. 0 .and. .not. overflow) then

C Create pointer from global to local orbitals
        do mu = 1,nbasis
          call GlobalToLocalOrb(mu,Node,Nodes,iptrgl(mu))
        enddo

        index = 0
        do ia = 1,natoms

C  determine how many LWF's depending on the atomic species 
C  (or the number of electrons)
          nelectr = qa (ia) + tiny
          if (ioptlwf .eq. 1) then
            indexi = ( ( nelectr + 2 ) / 2 )
          else if (ioptlwf .eq. 2) then
            indexi = ( ( nelectr ) / 2 )
          endif
	
C loop over LWF's centered on atom ia
	  do indexb = 1,indexi
	    index = index + 1

C initialize stuff...
            snor = 0.0d0
 
            do mu = lasto(ia-1)+1, lasto(ia)
              mul = iptrgl(mu)
              iorb = mu - lasto(ia-1)
              do in = 1,numc(mu)
                if (listc(in,mu) .eq. index) then
                  call initguess(ia,iorb,isa(ia),
     .                           nelectr,iseed,cg,Node)
                  if (mul.gt.0) c(in,mul) = cg
                  snor = snor + cg**2
                endif
              enddo
            enddo

C  Normalize LWF's if iopt = 0  .............................................
C  (normalize to one if functions are expected to be occupied, 0.5 if half
C  occupied and 0.1 if empty)
            fact = 1.0d0
            if (ioptlwf .eq. 1) then
              if (2*(nelectr/2) .eq. nelectr .and. indexb .eq. indexi) 
     .                       fact=sqrt(0.1d0)
              if (2*(nelectr/2) .ne. nelectr .and. indexb .eq. indexi) 
     .                       fact=sqrt(0.5d0)
            endif
            do mu = lasto(ia-1)+1, lasto(ia)
              mul = iptrgl(mu)
              if (mul.gt.0) then
                do in = 1,numc(mu)
                  if (listc(in,mu) .eq. index) then
                    c(in,mul) = c(in,mul) * fact / sqrt(snor)
                  endif
                enddo
              endif
            enddo

          enddo
        enddo
      endif

      if (index .gt. nbasis) then
        if (Node.eq.0) then
        write(6,*) 'cspa: ERROR: Num of LWFs larger than  basis set ',
     .    'size'
        write(6,*) '      ERROR: Increase basis set, or use less LWFs'
        endif
        stop
      endif

      if (2*index+tiny .lt. enum) then
        if (Node.eq.0) then
        write(6,*) 'cspa: ERROR: Too few LWFs for number of electrons'
        write(6,*) '      ERROR: Check if system charge is too large'
        endif
        stop
      endif

      if ((ioptlwf .eq. 2) .and. (nbands .ne. nqtot/2)) then
        if (Node.eq.0) then
        write(6,*) 'cspa: ERROR: Number of LWFs incorrectly calculated'
        write(6,*) '      ERROR: Something went wrong in generating the'
        write(6,*) '      ERROR: LWFs for the Ordejon-Mauri functional'
        endif
        stop
      endif

#ifdef MPI
C  Globalise maximum values
      ibuf1(1)=nctmax
      ibuf1(2)=ncmax
      call MPI_AllReduce(ibuf1,ibuf2,2,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)
      nctmax=ibuf2(1)
      ncmax=ibuf2(2)
      overflow=(ncmax.gt.maxnc)
#endif

C Deallocate local memory
      call memory('D','I',size(ilist),'cspa')
      deallocate(ilist)
      call memory('D','I',size(iptrgl),'cspa')
      deallocate(iptrgl)
      call memory('D','I',size(indexloc),'cspa')
      deallocate(indexloc)
      call memory('D','I',size(jan),'cspa')
      deallocate(jan)
      call memory('D','D',size(r2ij),'cspa')
      deallocate(r2ij)
      call memory('D','D',size(xij),'cspa')
      deallocate(xij)

      call timer('cspa',2)

      return
      end

      subroutine initguess(ia,iorb,is,ne,iseed,cg,Node)
C *****************************************************************************
C Routine to assign an initial guess for an atomic orbital iorb in a localized
C wave function centered on atom ia. 
C Assigns a random guess if the orbital belongs to the first 'zeta' of the
C atom (lowest energy shell of its angular momentum), and if the angular
C momentum is populated in the free atom. Otherwise, sets coefficient to cero.
C
C Written by P.Ordejon, November'96 
C lmax, lmaxs and nzls erased from the input by DSP, Aug. 1998.
C ******************************* INPUT ***************************************
C integer ia                   : Atom to which orbital belongs
C integer mu                   : Orbital index within atom ia
C integer is                   : Species index of atom ia
C integer ne                   : Number of electrons of atom ia
C integer iseed                : Seed for random number generator
C integer Node                 : Local node number
C ***************************** OUTPUT ****************************************
C real*8 cg                    : Initial guess for WF coefficient
C ***************************************************************************** 
C The following functions must exist:
C
C INTEGER FUNCTION LOMAXFIS(IS)
C    Returns the maximum angular momentum of orbitals
C Input:
C     INTEGER IS : Species index
C
C INTEGER FUNCTION NZTFL(IS,L)
C    Returns the number of different basis functions with the
C    same angular momentum L.
C Input:
C     INTEGER IS : Species index
C     INTEGER L  : Angular momentum
C
C ***************************************************************************** 
      use atmfuncs, only: lomaxfis, nztfl
      implicit none

      integer 
     .  ia,iorb,is,iseed,ne,Node


      double precision
     .  cg,ran3

      external ran3

C Internal variables .........................................................
      integer
     . index,iz,l,lmaxp,m
      
C .......................

C Initialize cg to cero
      cg = 0.0

C Find out angular momentum of orbital iorb

      index = 0
      do l = 0,lomaxfis(is)
        do iz = 1,nztfl(is,l)
          do m = -l,l
            index = index + 1
          enddo
        if (index .ge. iorb) goto 10
        enddo
      enddo
      if (Node.eq.0) then
        write(6,*) 'cspa: Error in orbital indexing in initguess'
      endif
      stop
10    continue

C Return if orbital is not the first zeta of its shell
      if (iz .ne. 1) return

C Assign initial guess.
C If 2 or less electrons, populate lowest s orbital
C If 8 or less electrons, populate lowest s and p  orbitals
C If 18 or less electrons, populate lowest s, p and d orbitals
C If 32 or less electrons, populate lowest s, p, d and f orbitals

      lmaxp = 0
      if (ne .le. 32) lmaxp = 3
      if (ne .le. 18) lmaxp = 2
      if (ne .le. 8) lmaxp = 1
      if (ne .le. 2) lmaxp = 0
      if (ne .gt. 32) then
        if (Node.eq.0) then
        write(6,*) 'cspa: Cannot build initial guess in initguess.'
        write(6,*) '      Reason: Too many electrons for this routine'
        endif
        stop
      endif 

      if (lmaxp .gt. lomaxfis(is)) then
        if (Node.eq.0) then
        write(6,*) 'cspa: Cannot build initial guess in initguess.'
        write(6,*) '      Reason: Max. angular moment for atom ',ia,
     .             '      is not large enough'
        endif
        stop
      endif 

      if (ne .gt. 32) then
        if (Node.eq.0) then
        write(6,*) 'cspa: Cannot build initial guess in initguess.'
        write(6,*) '      Too many valence electrons in atom ',ia
        endif
        stop
      endif 

      if (l .le. lmaxp) then
20      cg = (ran3(iseed) - 0.5d0) * 2.0d0
        if (abs(cg) .lt. 1d-5) goto 20
      endif

      return
      end
