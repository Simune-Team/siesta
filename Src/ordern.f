C $Id: ordern.f,v 1.16 1999/03/11 19:26:35 ordejon Exp $

      subroutine ordern(usesavelwf,ioptlwf,natoms,nbasis,lasto,
     .                  isa,qa,rcoor,rh,cell,xa,iscf,istep,itmax,
     .                  ftol,eta,enum,nhmax,numh,listh,h,s,
     .                  chebef,noeta,rcoorcp,beta,ipcheb,
     .                  dm,edm,Ecorrec)
C ****************************************************************************
C Order-N solver of the Hamiltonian problem.
C It calls the appropriate routines, and returns the
C density matrix to the calling program
C Uses the funcional of Kim et al. (PRB 52, 1640 (95))
C
C Written by P.Ordejon, October'96
C ****************************** INPUT *************************************
C logical usesavelwf           : True = try to use saved lwf from disk
C integer ioptlwf              : Build LWF's according to:
C                                0 = Read blindly from disk
C                                1 = Functional of Kim et al.
C                                2 = Functional of Ordejon-Mauri
C integer natoms               : Number of atoms
C integer nbasis               : Number of basis orbitals
C integer lasto(0:natoms)      : Index of last orbital of each atom
C integer isa(natoms)          : Species index of each atom
C real*8 qa(natoms)            : Neutral atom charge
C real*8 rcoor                 : Cutoff radius of Localized Wave Functions
C real*8 rh                    : Maximum cutoff radius of Hamiltonian matrix
C real*8 cell(3,3)             : Supercell vectors
C real*8 xa(3,natoms)          : Atomic coordinates
C integer iscf                 : SCF Iteration cycle (used to find 
C                                control vectors only if iscf=1)
C integer istep                : MD step
C integer itmax                : Maximum number of CG iterations
C real*8 ftol                  : Relative tolerance in CG minimization
C                                  (recomended: 1e-8)
C real*8 eta                   : Fermi level parameter of Kim et al.
C real*8 enum                   : Total number of electrons
C integer nhmax                : First dimension of listh and H, and maximum
C                                number of nonzero elements of each row of H
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of nonzero elements of each row of H)
C integer listh(nhmax,nbasis)  : Control vector of H matrix
C                                (list of nonzero elements of each row of H)
C real*8 h(nhmax,nbasis)       : Hamiltonian matrix (sparse)
C real*8 s(nhmax,nbasis)       : Overlap matrix (sparse)
C logical chebef               : Compute the chemical potential
C logical noeta                : Use computed Chem.pot. instead of eta
C real*8 rcoorcp               : Cutoff radius to compute Fermi level by
C                                projection.
C integer ipcheb               : Order of Chebishev expansion to compute Ef
C real*8 beta                  : Inverse Temperature for Chebishev expansion
C **************************** OUTPUT **************************************
C real*8 dm(nhmax,nbasis)      : Density matrix (sparse)
C real*8 edm(nhmax,nbasis)     : Energy Density matrix (sparse)
C real*8 Ecorrec               : Energy correction of Kim functional:
C                                eta * (etot-qs) , where qs is the charge
C                                of the Order-N solution
C **************************** BEHAVIOUR ***********************************
C If istep=1 and iscf=1, an initial guess is build. Otherwise, the LWF's of
C the former time steps are used,  extrapolated between MD steps.
C **************************************************************************** 

      implicit none

      include 'ordern.h'

      integer 
     .  ioptlwf, ipcheb, iscf, istep, itmax, natoms, 
     .  nbasis, nhmax

      integer
     .  isa(natoms), lasto(0:natoms), listh(nhmax,nbasis), 
     .  numh(nbasis)

      double precision
     .  beta, cell(3,3), dm(nhmax,nbasis), Ecorrec, edm(nhmax,nbasis), 
     .  enum, eta, ftol, h(nhmax,nbasis), qa(natoms), rcoor, rcoorcp,
     .  rh, s(nhmax,nbasis), xa(3,natoms)

      logical
     .  chebef, noeta, usesavelwf

C Internal variables ..........................................................
C maxoc = maximum number of atomic orbitals, needed to open numc, to be able
C         to find correct values for maxo and other dimensions
      integer maxoc
      parameter (maxoc = 20000)
      integer
     .  in, io, iopt, iord, iter, iterm, ncmax, nctmax, 
     .  nfmax, nftmax, nhijmax, nbands, nspin

      integer
     .  listc(maxnc,maxo), numc(maxoc),
     .  listcold(maxnc,maxo), numcold(maxoc)

      double precision
     .  aux(2,maxo), c(maxnc,maxo), cold(maxnc,maxo), 
     .  fe, qtot, xi(maxnc,maxo)

      real*4
     .  g(maxnc,maxo), hg(maxnc,maxo)

      double precision
     .  chpot,emax,emin,eV
C    .  gap,homo,lumo

      logical
     .  overflow, itest, found, frstme

      save 
     .  c, cold, frstme, iterm, listc, listcold, nbands, numc, numcold
     
      data frstme /.true./
      data iterm  / 0 /
C ...................

*     call timer( 'ordern', 1 )

      eV = 1.d0 / 13.60580d0

      overflow = .false.
C Check dimensions ......................................................
      call chkdim('ordern','maxoc',maxoc,nbasis,1) 
C ...................

C Print size of arrays ..................................................
      if (frstme) then
         call prmem( 0, 'ordern', 'aux',      'd', 2*maxo     )
         call prmem( 0, 'ordern', 'c',        'd', maxnc*maxo )
         call prmem( 0, 'ordern', 'cold',     'd', maxnc*maxo )
         call prmem( 0, 'ordern', 'g',        'd', maxnc*maxo )
         call prmem( 0, 'ordern', 'hg',       'd', maxnc*maxo )
         call prmem( 0, 'ordern', 'listc',    'i', maxnc*maxo )
         call prmem( 0, 'ordern', 'listcold', 'i', maxnc*maxo )
         call prmem( 0, 'ordern', 'xi',       'd', maxnc*maxo )
         call prmem( 0, 'ordern', ' ',        ' ', 0          )
      endif
C ............................

      write(6,"(/a,f12.4)") 'ordern: enum =', enum

C  Check if options are compatible
      if (ioptlwf .eq. 0 .and. (.not. usesavelwf)) then
        write(6,"(/a)") 'ordern: ERROR: You must use LWF files.'
        write(6,"(a)") '        If you choose ON.functional = Files'
        write(6,"(a)") '        Please set ON.UseSaveLWF = True'
        stop
      endif

C  If iscf = 1 (that is, if we are in a new MD step), find out initial
C  structure of localized wave functions, and initial guess .............
      if (iscf .eq. 1) then
        iopt = 1
        if (istep .eq. 1) then 
          if (usesavelwf) then
            call iolwf( 'read', maxnc, maxo, nbasis, 1,
     .                  numcold, listcold, c, found )
            if (found) then
C             Find out number of bands
              nbands = 0
              do io = 1,nbasis
                do in = 1,numcold(io)
                  nbands = max( nbands, listcold(in,io) )
                enddo
              enddo
            else
              iopt = 0
*             write(6,"(/a)")'ordern: You must have LWF files available'
*             write(6,"(a)")'        if you choose ON.functional = Files'
*             stop
            endif
          else
            iopt = 0
          endif
        endif
        if (nbasis .gt. maxo)  overflow = .true.

        if (ioptlwf .ne. 0) then
C   This call was modified by DSP, Aug. 1998.
          call cspa(ioptlwf,iopt,natoms,nbasis,
     .        lasto,isa,
     .        qa,rcoor,rh,cell,xa,nhmax,numh,listh,maxnc,
     .        c,numc,listc,ncmax,nctmax,nfmax,nftmax,nhijmax,nbands,
     .        overflow)

C         Check that dimensios are large enough 
          if (nbasis   .gt. maxo)    overflow = .true.
          if (nhmax    .gt. maxnh)   overflow = .true.
          if (nbands   .gt. maxlwf)  overflow = .true.
          if (ncmax    .gt. maxnc)   overflow = .true.
          if (nctmax   .gt. maxnct)  overflow = .true.
          if (nfmax    .gt. maxnf)   overflow = .true.
          if (nftmax   .gt. maxnft)  overflow = .true.
          if (nhijmax  .gt. maxnhij) overflow = .true.
          if (max(nfmax,nhijmax) .gt. maxnhf) overflow = .true.
   
C         Print good dimensions
          if (overflow) then
          open( 1, file='ordern.h', status='unknown' )
          write(1,'(A)') 'C ordern: Dimension parameters for ordern'
          write(1,'(6X,A)')'integer maxnc,maxnct,maxnf,maxnft,maxnh'
          write(1,'(6X,A)')'integer maxnhij,maxnhf,maxlwf,maxo'
          write(1,'(6X,A,I10,A)') 'parameter ( maxo   = ',  nbasis, ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxnh  = ',  nhmax,  ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxlwf = ',  nbands, ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxnc  = ',   ncmax, ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxnct = ',  nctmax, ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxnf  = ',   nfmax, ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxnft = ',  nftmax, ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxnhij= ', nhijmax, ' )'
          write(1,'(6X,A,I10,A)') 'parameter ( maxnhf = ',  
     .                                         max(nfmax,nhijmax),  ' )'
          close(1)
          write(6,'(/,a,/)') 'ordern: BAD DIMENSIONS. RECOMPILE.'
          stop 'ordern: BAD DIMENSIONS. RECOMPILE.'
          endif

        endif
      endif
C ..................

C Calculate Chemical Potential, Max and Min eigenvalues, energy gap
C and HOMO and LUMO levels ..........

      if (chebef) then

        call timer( 'chempot', 1 )
        call chempot(h,s,listh,numh,rcoorcp,ipcheb,beta,lasto,
     .               cell,xa,enum,nbasis,natoms,nhmax,
     .               chpot,emax,emin)
C     . ,gap,homo,lumo)
        call timer( 'chempot', 2 )

        eV  = 1.d0 / 13.60580d0

        write(6,'(a,f8.4,a)') 'ordern:   Chemical Potential = ',
     .                         chpot/eV,' eV'
        write(6,'(a,f8.4,a)') 'ordern:   Maximum Eigenvalue = ',
     .                         emax/eV,' eV'
        write(6,'(a,f8.4,a)') 'ordern:   Minimum Eigenvalue = ',
     .                         emin/eV,' eV'
c      write(6,*) 'Gap energy         = ',gap/eV,' eV'
c      write(6,*) 'HOMO               = ',homo/eV,' eV'
c      write(6,*) 'LUMO               = ',lumo/eV,' eV'

        if (noeta) eta=chpot
      endif


C  Extrapolate wave funcions from those of former time steps if 
C  iscf = 1 (that is, if we are in a new MD step) ...................
      if (iscf .eq. 1)  then
        write(6,"(/a,i3)") 'ordern: ioptlwf =',ioptlwf
        if (ioptlwf .eq. 0) then
          do io = 1,nbasis
            numc(io) = numcold(io)
            do in = 1,numcold(io)
              listc(in,io) = listcold(in,io)
*              c(in,io) = cold(in,io)
            enddo
          enddo
        endif
        iord = 1
	if (iterm .gt. 50) iord = 0
C       If LWF's have just been read from disk, 
C       call extrapol with istep = 2 and iord = 1
C       to make it update the structure of c, if needed
        itest = .false.
        if (istep.eq.1 .and. usesavelwf .and. found) then
          istep = 2
          iord = 0
          itest = .true.
        endif
        nspin = 1
        call extrapol(istep,iord,nspin,nbasis,maxo,maxnc,numc,listc,
     .                aux,numcold,listcold,cold,c)
C       If LWF's have just been read, restore istep
        if (itest) istep = 1
        itest = .false.
      endif
C .................

C Call the CG routines ...............................................
      if (iscf .eq. 1) iterm = 0
      call cgwf(iscf,itmax,ftol,eta,enum,nbasis,nbands,
     .           nhmax,numh,listh,maxnc,numc,listc,h,s,c,
     .           g,hg,xi,
     .           fe,iter,dm,edm)
      iterm = max(iterm, iter)
C ..................

C Calculate correction to the total energy from Kim's functional .....
C First calculate total charge of the solution
      qtot = 0.d0
      do io = 1,nbasis
        do in = 1,numh(io)
          qtot = qtot + dm(in,io) * s(in,io)
        enddo
      enddo

      Ecorrec = eta * (enum - qtot)
      write(6,"(a,f12.4)") 'ordern: qtot (after  DM normalization) = ',
     .     qtot
C ..................

C Save LWF's to disk .................................................
      call iolwf( 'write', maxnc, maxo, nbasis, 1,
     .             numc, listc, c, found )
C ....................

      frstme = .false.
*     call timer( 'ordern', 2 )
      return
      end


