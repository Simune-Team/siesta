      subroutine ordern(usesavelwf,ioptlwf,natoms,nbasis,lasto,
     .                  isa,qa,rcoor,rh,cell,xa,iscf,istep,itmax,
     .                  ftol,eta,enum,maxuo,nhmax,numh,listh,h,s,
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
C integer maxuo                : Second dimension of H, S, dm and edm.
C integer nhmax                : First dimension of listh and H, and maximum
C                                number of nonzero elements of each row of H
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of nonzero elements of each row of H)
C integer listh(nhmax,nbasis)  : Control vector of H matrix
C                                (list of nonzero elements of each row of H)
C real*8 h(nhmax,maxuo)        : Hamiltonian matrix (sparse)
C real*8 s(nhmax,maxuo)        : Overlap matrix (sparse)
C logical chebef               : Compute the chemical potential
C logical noeta                : Use computed Chem.pot. instead of eta
C real*8 rcoorcp               : Cutoff radius to compute Fermi level by
C                                projection.
C integer ipcheb               : Order of Chebishev expansion to compute Ef
C real*8 beta                  : Inverse Temperature for Chebishev expansion
C **************************** OUTPUT **************************************
C real*8 dm(nhmax,maxuo)       : Density matrix (sparse)
C real*8 edm(nhmax,maxuo)      : Energy Density matrix (sparse)
C real*8 Ecorrec               : Energy correction of Kim functional:
C                                eta * (etot-qs) , where qs is the charge
C                                of the Order-N solution
C **************************** BEHAVIOUR ***********************************
C If istep=1 and iscf=1, an initial guess is build. Otherwise, the LWF's of
C the former time steps are used,  extrapolated between MD steps.
C **************************************************************************** 

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
     .  ioptlwf, ipcheb, iscf, istep, itmax, natoms, 
     .  nbasis, nhmax, maxuo

      integer
     .  isa(natoms), lasto(0:natoms), listh(nhmax,nbasis), 
     .  numh(nbasis)

      double precision
     .  cell(3,3), dm(nhmax,maxuo), Ecorrec, edm(nhmax,maxuo), 
     .  enum, eta, ftol, h(nhmax,maxuo), qa(natoms), rcoor, rh, 
     .  s(nhmax,maxuo), xa(3,natoms), beta, rcoorcp

      logical
     .  chebef, noeta, usesavelwf

C Internal variables ..........................................................
      integer, save ::
     .  maxo, maxoloc, maxlwf, maxnc, maxnct, maxnf, maxnft,
     .  maxnhij

      integer
     .  in, io, iio, iopt, iord, iter, iterm, ncmax, nctmax, 
     .  nfmax, nftmax, nhijmax, nbands, nspin, Node, Nodes,
     .  nbasisloc, nbasismax, nbandsloc

      integer, dimension(:,:), allocatable, save ::
     .  listc, listcold

      integer, dimension(:), allocatable, save ::
     .  numc, numcold

      double precision
     .  fe, qtot, chpot,emax,emin,eV

      real*4, dimension(:,:), allocatable, save :: 
     .  g, hg

      double precision, dimension(:,:), allocatable, save ::
     .  xi

      double precision, dimension(:,:), allocatable, save ::
     .  c, cold

      logical
     .  overflow, itest, found 

#ifdef MPI
      integer
     .  MPIerror
      double precision
     .  rtmp
#endif

      save 
     .  iterm, nbands
     
      data iterm   / 0 /
      data maxo    / 1 /
      data maxoloc / 1 /
      data maxlwf  / 1 /
      data maxnc   / 1 /
      data maxnct  / 1 /
      data maxnf   / 1 /
      data maxnft  / 1 /
      data maxnhij / 1 /
C ...................

C Get Node number and local number of orbitals in unit cell
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
      call GetNodeOrbs(nbasis,Node,Nodes,nbasisloc)
      call MPI_AllReduce(nbasisloc,nbasismax,1,MPI_integer,
     .  MPI_Max,MPI_Comm_World,MPIerror)
#else
      Node = 0
      Nodes = 1
      nbasisloc = nbasis
      nbasismax = nbasis
#endif

*     call timer( 'ordern', 1 )

      eV = 1.d0 / 13.60580d0

      overflow = .false.
C ...................

C Allocate local memory - never deallocated as contents must be saved
      maxo = nbasis
      maxoloc = nbasismax
      if (.not.allocated(numc)) then
        allocate(numc(maxo))
        call memory('A','I',maxo,'ordern')
      endif
      if (.not.allocated(numcold)) then
        allocate(numcold(maxo))
        call memory('A','I',maxo,'ordern')
      endif
      if (.not.allocated(listc)) then
        allocate(listc(maxnc,maxo))
        call memory('A','I',maxnc*maxo,'ordern')
      endif
      if (.not.allocated(listcold)) then
        allocate(listcold(maxnc,maxo))
        call memory('A','I',maxnc*maxo,'ordern')
      endif
      if (.not.allocated(c)) then
        allocate(c(maxnc,maxoloc))
        call memory('A','D',maxnc*maxoloc,'ordern')
      endif
      if (.not.allocated(cold)) then
        allocate(cold(maxnc,maxoloc))
        call memory('A','D',maxnc*maxoloc,'ordern')
      endif

      if (Node.eq.0) then
        write(6,"(/a,f12.4)") 'ordern: enum =', enum
      endif

C  Check if options are compatible
      if (ioptlwf .eq. 0 .and. (.not. usesavelwf)) then
        if (Node.eq.0) then
          write(6,"(/a)") 'ordern: ERROR: You must use LWF files.'
          write(6,"(a)") '        If you choose ON.functional = Files'
          write(6,"(a)") '        Please set ON.UseSaveLWF = True'
        endif
        stop
      endif

C  If iscf = 1 (that is, if we are in a new MD step), find out initial
C  structure of localized wave functions, and initial guess .............
      if (iscf .eq. 1) then
   10   iopt = 1
        if (istep .eq. 1) then 
          if (usesavelwf) then
            call iolwf( 'read', maxnc, maxo, maxoloc, nbasis, 1,
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
            endif
          else
            iopt = 0
          endif
        endif
        if (nbasis .gt. maxo)  overflow = .true.
#ifdef MPI
        call GetNodeOrbs(nbasis,Node,Nodes,nbasisloc)
        call MPI_AllReduce(nbasisloc,nbasismax,1,MPI_integer,
     .    MPI_Max,MPI_Comm_World,MPIerror)
#else
        nbasisloc = nbasis
        nbasismax = nbasis
#endif
        if (nbasisloc .gt. maxoloc)  overflow = .true.

C  If overflow re-dimension arrays and repeat previous step
        if (overflow) then
          call memory('D','I',size(numc),'ordern')
          deallocate(numc)
          call memory('D','I',size(numcold),'ordern')
          deallocate(numcold)
          call memory('D','I',size(listc),'ordern')
          deallocate(listc)
          call memory('D','I',size(listcold),'ordern')
          deallocate(listcold)
          call memory('D','D',size(c),'ordern')
          deallocate(c)
          call memory('D','D',size(cold),'ordern')
          deallocate(cold)
          maxo = nbasis
          maxoloc = nbasismax
          allocate(numc(maxo))
          call memory('A','I',maxo,'ordern')
          allocate(numcold(maxo))
          call memory('A','I',maxo,'ordern')
          allocate(listc(maxnc,maxo))
          call memory('A','I',maxnc*maxo,'ordern')
          allocate(listcold(maxnc,maxo))
          call memory('A','I',maxnc*maxo,'ordern')
          allocate(c(maxnc,maxoloc))
          call memory('A','D',maxnc*maxoloc,'ordern')
          allocate(cold(maxnc,maxoloc))
          call memory('A','D',maxnc*maxoloc,'ordern')
          overflow = .false.
          goto 10
        endif

        nbandsloc = ((nbands-1)/Nodes) + 1

        if (ioptlwf .ne. 0) then
C   This call was modified by DSP, Aug. 1998.
   20     call cspa(ioptlwf,iopt,natoms,nbasis,nbasisloc,lasto,isa,
     .        qa,enum,rcoor,rh,cell,xa,nhmax,numh,listh,maxnc,
     .        c,numc,listc,ncmax,nctmax,nfmax,nftmax,nhijmax,nbands,
     .        overflow)

          nbandsloc = ((nbands-1)/Nodes) + 1

C         Check that dimensions are large enough 
          if (nbasis   .gt. maxo)    overflow = .true.
          if (nbasismax.gt. maxoloc) overflow = .true.
          if (nbandsloc.gt. maxlwf)  overflow = .true.
          if (ncmax    .gt. maxnc)   overflow = .true.
   
C         Print good dimensions
          if (overflow) then
            if (nbasis .gt. maxo) then
              call memory('D','I',size(numc),'ordern')
              deallocate(numc)
              call memory('D','I',size(numcold),'ordern')
              deallocate(numcold)
              allocate(numc(nbasis))
              call memory('A','I',nbasis,'ordern')
              allocate(numcold(nbasis))
              call memory('A','I',nbasis,'ordern')
            endif
            if (nbasis .gt. maxo .or. ncmax .gt. maxnc) then
              call memory('D','I',size(listc),'ordern')
              deallocate(listc)
              call memory('D','I',size(listcold),'ordern')
              deallocate(listcold)
              allocate(listc(ncmax,nbasis))
              call memory('A','I',ncmax*nbasis,'ordern')
              allocate(listcold(ncmax,nbasis))
              call memory('A','I',ncmax*nbasis,'ordern')
            endif
            if (nbasismax .gt. maxoloc .or. ncmax .gt. maxnc) then
              call memory('D','D',size(c),'ordern')
              deallocate(c)
              call memory('D','D',size(cold),'ordern')
              deallocate(cold)
              allocate(c(ncmax,nbasismax))
              call memory('A','D',ncmax*nbasismax,'ordern')
              allocate(cold(ncmax,nbasismax))
              call memory('A','D',ncmax*nbasismax,'ordern')
            endif
            maxo = nbasis
            maxoloc = nbasismax
            maxlwf = nbandsloc
            maxnc = ncmax
            overflow = .false.
            goto 20
          endif
          maxnf = nfmax
          maxnft = nftmax
          maxnhij = nhijmax
          maxnct = nctmax

        endif
      endif
C ..................

C Calculate Chemical Potential, Max and Min eigenvalues, energy gap
C and HOMO and LUMO levels ..........

      if (chebef) then

        call timer( 'chempot', 1 )
        call chempot(h,s,listh,numh,rcoorcp,ipcheb,beta,lasto,
     .               cell,xa,enum,nbasis,nbasisloc,natoms,nhmax,
     .               chpot,emax,emin)
C     . ,gap,homo,lumo)
        call timer( 'chempot', 2 )

        if (Node.eq.0) then
        write(6,'(a,f8.4,a)') 'ordern:   Chemical Potential = ',
     .                         chpot/eV,' eV'
        write(6,'(a,f8.4,a)') 'ordern:   Maximum Eigenvalue = ',
     .                         emax/eV,' eV'
        write(6,'(a,f8.4,a)') 'ordern:   Minimum Eigenvalue = ',
     .                         emin/eV,' eV'
c       write(6,*) 'Gap energy         = ',gap/eV,' eV'
c       write(6,*) 'HOMO               = ',homo/eV,' eV'
c       write(6,*) 'LUMO               = ',lumo/eV,' eV'
        endif

        if (noeta) eta=chpot
      endif

C  Extrapolate wave funcions from those of former time steps if 
C  iscf = 1 (that is, if we are in a new MD step) ...................
      if (iscf .eq. 1)  then
        if (Node.eq.0) then
          write(6,"(/a,i3)") 'ordern: ioptlwf =',ioptlwf
        endif
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
        call extrapol(istep,iord,nspin,nbasisloc,maxo,maxoloc,
     .                maxnc,numc,listc,numcold,listcold,
     .                cold,c)
C       If LWF's have just been read, restore istep
        if (itest) istep = 1
        itest = .false.
      endif
C .................

C Call the CG routines ...............................................
      if (iscf .eq. 1) iterm = 0
      allocate(xi(maxnc,maxoloc))
      call memory('A','D',maxnc*maxoloc,'ordern')
      allocate(hg(maxnc,maxoloc))
      call memory('A','S',maxnc*maxoloc,'ordern')
      allocate(g(maxnc,maxoloc))
      call memory('A','S',maxnc*maxoloc,'ordern')
      call cgwf(iscf,itmax,ftol,eta,enum,nbasis,nbands,
     .           nhmax,numh,listh,maxnc,numc,listc,h,s,c,
     .           g,hg,xi,fe,iter,dm,edm,nbasisloc,
     .           maxnf,maxnft,maxnhij,maxnct)
      call memory('D','D',size(xi),'ordern')
      deallocate(xi)
      call memory('D','S',size(hg),'ordern')
      deallocate(hg)
      call memory('D','S',size(g),'ordern')
      deallocate(g)
      iterm = max(iterm, iter)
C ..................

C Calculate correction to the total energy from Kim's functional .....
C First calculate total charge of the solution
      qtot = 0.d0
      do io = 1,nbasisloc
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        do in = 1,numh(iio)
          qtot = qtot + dm(in,io) * s(in,io)
        enddo
      enddo

#ifdef MPI
C Global reduction of qtot
      call MPI_AllReduce(qtot,rtmp,1,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      qtot = rtmp
#endif

      Ecorrec = eta * (enum - qtot)
      if (Node.eq.0) then
      write(6,"(a,f12.4)") 'ordern: qtot (after  DM normalization) = ',
     .   qtot
      endif
C ..................

C Save LWF's to disk .................................................
      call iolwf( 'write', maxnc, maxo, maxoloc, nbasis, 1,
     .             numc, listc, c, found )
C ....................

*     call timer( 'ordern', 2 )
      return
      end


