C $Id: atomlwf.f,v 1.10.2.2 1999/06/08 10:40:01 ordejon Exp $

      subroutine cspa(ioptlwf,iopt,natoms,nbasis,lasto,isa,
     .        qa,enum,rcoor,rh,cell,xa,nhmax,numh,listh,maxnc,
     .        c,numc,listc,ncmax,nctmax,nfmax,nftmax,nhijmax,nbands,
     .        overflow)
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
C integer nfmax             : Maximum number of nonzero elements 
C                             of each row of F
C integer nftmax            : Maximum number of nonzero elements 
C                             of eaxh column of F
C integer nhijmax           : Maximum number of nonzero elements 
C                             of each row of Hij
C integer nbands            : Number of LWF's
C ************************** INPUT AND OUTPUT ********************************
C logical overflow          : True: Maximum dimensions too small
C ****************************************************************************

      implicit none

      integer 
     .  iopt,ioptlwf,maxnc,natoms,nbands,nbasis,ncmax,nctmax,
     .  nfmax,nftmax,nhmax,nhijmax

      integer 
     .  isa(natoms),lasto(0:natoms),
     .  numc(nbasis),listc(maxnc,nbasis),numh(nbasis),
     .  listh(nhmax,nbasis)

      double precision
     .  cell(3,3),enum,qa(natoms),rcoor,rh,xa(3,natoms),
     .  c(maxnc,nbasis)
 
      logical
     .  overflow

C  Internal variables .......................................................
C  nnmax = maximum number of neighbour atoms within rcoor
      integer nnmax
      parameter (nnmax = 1000)

C  nbmax = maximum number of basis orbitals
      integer nbmax
      parameter (nbmax = 10000)

C  namax = maximum number of atoms
      integer namax
      parameter (namax = 2000)

      integer
     . alist(namax),i,ia,ilist(nbmax),imu,
     . in,index,indexa,indexb,indexi,indexj,indexloc(nnmax),iorb,iseed,
     . j,ja,jan(nnmax),jj,k,mu,nct,nelectr,nf,nm,nna,norb,nqtot,nu,
     . numft(nbmax),numloc,nhij

      double precision
     . cg,fact,qtot,r,r2ij(nnmax),rmax,rr(3),rrmod,
     . snor,tiny,xij(3,nnmax)
C ...........................................................................

      tiny = 1.d-10

C Check that iopt is correct ................................................
      if (iopt .ne. 0 .and. iopt .ne. 1) then
        write(6,*) 'cspa: Wrong iopt in cspa'
        stop
      endif
C ..................

C Check dimensions ..........................................................
      if (natoms .gt. namax) then
        write(6,*) 'cspa: Wrong namax; Must be at least ',natoms
        stop
      endif
      if (nbasis .gt. nbmax) then
        write(6,*) 'cspa: Wrong nbmax; Must be at least ',nbasis
        stop
      endif
C ..................

C Initialize some stuff ......................................................

      do ia = 1,natoms
        alist(ia) = 0
      enddo

      do mu = 1,nbasis
        numc(mu) = 0
        numft(mu) = 0
      enddo

      if (iopt .eq. 0) then
        do mu = 1,nbasis
          do i = 1,maxnc
            c(i,mu) = 0.0
          enddo
        enddo
      endif

      iseed = 17

      ncmax   = 0
      nctmax  = 0
      nfmax   = 0
      nftmax  = 0
      nhijmax = 0
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
      nqtot = enum + tiny
      if (abs(nqtot-enum+tiny) .gt. 1e-3) then
        write(6,*) 'cspa: WARNING: total charge non integer:',enum
      endif
      if (2*(nqtot/2) .ne. nqtot) then
        write(6,*) 'cspa: ERROR: Wrong total charge; odd charge:',qtot
        write(6,*) '      ERROR: Charge must be EVEN to use Order-N ',
     .    'option'
        stop
      endif
C ..................

C Check if system is charged..............
      qtot = 0.d0
      do ia = 1,natoms
        qtot = qtot + qa(ia)
      enddo
      if (abs(qtot-enum) .gt. 1.0d-4) then
        if (ioptlwf .eq. 2) then
          write(6,*) 'cspa: ERROR: Charged systems can not be done'
          write(6,*) '      ERROR: with the Ordejon-Mauri functional'
          stop
        endif
      endif
C ..................

C Build control vectors of sparse LWF 
C loop over the localized wave funcions (centered at atoms)..................

C Initialize routine for neighbour search
        nna = nnmax
        call neighb(cell,rcoor,natoms,xa,0,0,nna,jan,xij,r2ij)

      index = 0
      do ia = 1,natoms

        if (2.*rcoor .lt. rmax) then
C Look for neighbours of atom ia
          nna = nnmax
          call neighb(cell,rcoor,natoms,xa,ia,0,nna,jan,xij,r2ij)
          if (nna .gt. nnmax) then
            write(6,*) 'cspa: nnmax = ',nnmax,' ; should be at least '
     .                  ,nna
            stop
          endif
        else
          nna = natoms
          do  jj = 1,natoms
            jan(jj) = jj
          enddo
        endif


C  determine how many LWF's depending on the atomic species 
C  (or the number of electrons)
        nelectr = qa (ia) + tiny
        if (abs(nelectr - qa(ia) + tiny) .gt. 1e-3) then
          write(6,*) 'cspa: Wrong atomic charge for atom ',ia
          write(6,*) '      qa = ',qa(ia),'  Must be an integer!!'
          stop
        endif
        if (ioptlwf .eq. 1) then
          indexi = ( ( nelectr + 2 ) / 2 )
        else if (ioptlwf .eq. 2) then
          if ( (nelectr/2)*2 .ne. nelectr) then
            write(6,*) 'cspa: Wrong Order-N functional option in cspa.'
            write(6,*) '      You can only use the functional of'
            write(6,*) '      Ordejon-Mauri for atoms with an even'
            write(6,*) '      number of electrons.'
            stop
          endif
          indexi = ( ( nelectr ) / 2 )
        else
          write(6,*) 'cspa: Wrong functional option in cspa'
          stop
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
          snor = 0.0
          nf = 0
          do nu = 1,nbasis
            ilist(nu) = 0
          enddo
C ...
 
c  loop over the neighbors of ia within rcoor
          do  jj = 1,nna
            ja = jan(jj)
            norb = lasto(ja)-lasto(ja-1)

c  check if ja has already been included in current lwf
            do indexa = 1,numloc
              if (ja .eq. indexloc(indexa)) goto 30
            enddo
            numloc = numloc + 1
            indexloc(numloc) = ja

C  loop over orbitals of ja
            do iorb = 1,norb
              mu = iorb + lasto(ja-1)
              nm = numc(mu)
              numc(mu) = nm + 1
              if (numc(mu) .gt. maxnc) overflow = .true.
              if (.not.overflow) listc(nm+1,mu) = index
              nct = nct + 1

C Find out structure of F and Ft matrices
C   find orbitals which interact with mu
              do imu = 1,numh(mu)
                nu = listh(imu,mu)
                if (ilist(nu) .eq. 0) then
                  ilist(nu) = 1
                  numft(nu) = numft(nu) + 1
                  nf = nf + 1
                endif
              enddo

C  Assign ramdom guess for orbitals in atom ia if iopt = 0
              if (.not.overflow) then
                if (iopt .eq. 0) then
                  if (ja .eq. ia) then
                    call initguess(ia,iorb,isa(ia),
     .                             nelectr,iseed,cg)
                    c(nm+1,mu) = cg
                    snor = snor + cg**2
                  endif
                endif
              endif
            enddo
30        enddo

C Normalize LWF's if iopt = 0  .............................................
C  (normalize to one if functions are expected to be occupied, 0.5 if half
C   occupied and 0.1 if empty)
          if (.not.overflow)  then
            if (iopt .eq. 0) then
              fact = 1.0
              if (ioptlwf .eq. 1) then
                if (2*(nelectr/2) .eq. nelectr .and. indexb .eq. indexi) 
     .                         fact=sqrt(0.1)
                if (2*(nelectr/2) .ne. nelectr .and. indexb .eq. indexi) 
     .                         fact=sqrt(0.5)
              endif
              do mu = lasto(ia-1)+1, lasto(ia)
                do in = 1,numc(mu)
                  if (listc(in,mu) .eq. index) then
                    c(in,mu) = c(in,mu) * fact / sqrt(snor)
                  endif
                enddo
              enddo
            endif
          endif
C ..........

          nctmax = max ( nctmax , nct )
          nfmax = max ( nfmax , nf )
          
        enddo
      enddo

      do mu = 1,nbasis
        nftmax = max ( nftmax , numft(mu) )
        ncmax  = max ( ncmax  , numc(mu) )
      enddo

      nbands = index

      if (index .gt. nbasis) then
        write(6,*) 'cspa: ERROR: Num of LWFs larger than  basis set ',
     .    'size'
        write(6,*) '      ERROR: Increase basis set, or use less LWFs'
        stop
      endif

      if (2*index+tiny .lt. enum) then
        write(6,*) 'cspa: ERROR: Too few LWFs for number of electrons'
        write(6,*) '      ERROR: Check if system charge is too large'
        stop
      endif

      if ((ioptlwf .eq. 2) .and. (nbands .ne. nqtot/2)) then
        write(6,*) 'cspa: ERROR: Number of LWFs incorrectly calculated'
        write(6,*) '      ERROR: Something went wrong in generating the'
        write(6,*) '      ERROR: LWFs for the Ordejon-Mauri functional'
        stop
      endif


C Find out sparse dimensions of Hij
C loop over the localized wave funcions (centered at atoms)..................

C Maximum interacion range between LWF's
      r = 2.0 * rcoor + rh

      if (2*r .ge. rmax) then
        nhijmax = nbands
      else
        nna = nnmax
        call neighb(cell,r,natoms,xa,0,0,nna,jan,xij,r2ij)

        index = 0
        do ia = 1,natoms
  
C Look for neighbours of atom ia within maximum interaction range
          nna = nnmax
          call neighb(cell,r,natoms,xa,ia,0,nna,jan,xij,r2ij)
          if (nna .gt. nnmax) then
            write(6,*) 'cspa: nnmax = ',nnmax,' ; should be at least '
     .                 ,nna
            stop
          endif

          nhij = 0
c  loop over the neighbors of ia within rcoor
          do  jj = 1,nna
            ja=jan(jj)
            alist(ja) = 0
          enddo
          do  jj = 1,nna
            ja = jan(jj)
            if (alist(ja) .eq. 1) goto 20
            alist(ja) = 1

C  determine how many LWF's centered in ja, depending on the atomic species 
C  (or the number of electrons)
            nelectr = qa (ja) + tiny
            if (ioptlwf .eq. 1) then
              indexj = ( ( nelectr + 2 ) / 2 )
            else if (ioptlwf .eq. 2) then
              indexj = ( ( nelectr / 2 ) )
            else
              write(6,*) 'cspa: Wrong functional option in cspa'
              stop
            endif

            nhij = nhij + indexj
20          continue
          enddo
          do  jj = 1,nna
            ja=jan(jj)
            alist(ja) = 0
          enddo
          nhijmax = max ( nhijmax , nhij )
        enddo
      endif
C ............

      return
      end



      subroutine initguess(ia,iorb,is,ne,iseed,cg)
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

      implicit none

      integer 
     .  ia,iorb,is,iseed,ne


      double precision
     .  cg,ran3

      external ran3

C Internal variables .........................................................
      integer
     . index,iz,l,lmaxp,m, lomaxfis,nztfl
      
      external 
     . lomaxfis, nztfl
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
      write(6,*) 'cspa: Error in orbital indexing in initguess'
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
        write(6,*) 'cspa: Cannot build initial guess in initguess.'
        write(6,*) '      Reason: Too many electrons for this routine'
        stop
      endif 

      if (lmaxp .gt. lomaxfis(is)) then
        write(6,*) 'cspa: Cannot build initial guess in initguess.'
        write(6,*) '      Reason: Max. angular moment for atom ',ia,
     .             '      is not large enough'
        stop
      endif 

      if (ne .gt. 32) then
        write(6,*) 'cspa: Cannot build initial guess in initguess.'
        write(6,*) '      Too many valence electrons in atom ',ia
        stop
      endif 

      if (l .le. lmaxp) then
20      cg = (ran3(iseed) - 0.5) * 2.0
        if (abs(cg) .lt. 1d-5) goto 20
      endif

      return
      end
