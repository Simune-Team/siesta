      module m_ordern
!
!     We need an explicit interface since assumed-shape arrays are used
!
      CONTAINS

      subroutine ordern(usesavelwf,ioptlwf,natoms,nbasis,lasto,
     .                  isa,qa,rcoor,rh,cell,xa,iscf,istep,itmax,
     .                  ftol,eta,enum,maxuo,maxnh,numh,listhptr,
     .                  listh,h,s,chebef,noeta,rcoorcp,beta,ipcheb,
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
C integer maxnh                : First dimension of listh and H, and maximum
C                                number of nonzero elements of H
C integer numh(nbasis)         : Control vector of H matrix
C                                (number of nonzero elements of each row of H)
C integer listh(maxnh)         : Control vector of H matrix
C                                (list of nonzero elements of each row of H)
C real*8 h(maxnh)              : Hamiltonian matrix (sparse)
C real*8 s(maxnh)              : Overlap matrix (sparse)
C logical chebef               : Compute the chemical potential
C logical noeta                : Use computed Chem.pot. instead of eta
C real*8 rcoorcp               : Cutoff radius to compute Fermi level by
C                                projection.
C integer ipcheb               : Order of Chebishev expansion to compute Ef
C real*8 beta                  : Inverse Temperature for Chebishev expansion
C **************************** OUTPUT **************************************
C real*8 dm(maxnh)             : Density matrix (sparse)
C real*8 edm(maxnh)            : Energy Density matrix (sparse)
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
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      double precision, intent(INOUT) :: eta   !!  ?????????
      integer, intent(inout)          :: istep

      integer, intent(in)           ::
     .  ioptlwf, ipcheb, iscf, itmax, natoms, 
     .  nbasis, maxnh, maxuo

      integer, intent(in)           ::
     .  isa(natoms), lasto(0:natoms), listh(maxnh), 
     .  numh(:), listhptr(:)

      double precision, intent(out) :: dm(maxnh), Ecorrec, edm(maxnh)
      double precision, intent(in)  ::
     .  cell(3,3), 
     .  enum, ftol, h(maxnh), qa(natoms), rcoor, rh, 
     .  s(maxnh), xa(3,natoms), beta, rcoorcp

      logical, intent(in)           :: chebef, noeta, usesavelwf

C Internal variables ..........................................................
      integer, save ::
     .  maxo, maxoloc, maxlwf, maxnc, maxnct

      integer
     .  in, io, iio, iopt, iord, iter, iterm, ncmax, nctmax, 
     .  nh, nbands, nspin, Node, Nodes,
     .  nbasisloc, nbasismax, nbandsloc, ntmp, BNode

      integer, dimension(:,:), allocatable, save ::
     .  listc, listcold, ntmp2

      integer, dimension(:), allocatable, save ::
     .  numc, numcold, ntmp1

      double precision
     .  fe, qtot, chpot,emax,emin,eV

      real*4, dimension(:,:), allocatable, save :: 
     .  g, hg

      double precision, dimension(:,:), allocatable, save ::
     .  xi

      double precision, dimension(:,:), allocatable, save ::
     .  c, cold, rtmp2

      logical
     .  overflow, itest, found, writelwf

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
        allocate(listcold(maxnc,maxoloc))
        call memory('A','I',maxnc*maxoloc,'ordern')
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
            call iolwf( 'read', maxnc, ncmax, maxoloc, nbasis, 1,
     .                  numcold, listcold, cold, found )
            if (found.and.ncmax.le.maxnc) then
C             Find out number of bands
              nbands = 0
              do io = 1,nbasis
                call GlobalToLocalOrb(io,Node,Nodes,iio)
                if (iio.gt.0) then
                  do in = 1,numcold(io)
                    nbands = max( nbands, listcold(in,iio) )
                  enddo
                endif
              enddo
C Copy Cold to C
              c = cold
            else
              iopt = 0
            endif
          else
            iopt = 0
            ncmax = maxnc
          endif
        else 
          ncmax = maxnc
        endif
        if (nbasis .gt. maxo)  overflow = .true.
#ifdef MPI
        call GetNodeOrbs(nbasis,Node,Nodes,nbasisloc)
        call MPI_AllReduce(nbasisloc,nbasismax,1,MPI_integer,
     .    MPI_Max,MPI_Comm_World,MPIerror)
        call MPI_AllReduce(nbands,ntmp,1,MPI_integer,
     .    MPI_Max,MPI_Comm_World,MPIerror)
        nbands = ntmp
        call MPI_AllReduce(ncmax,ntmp,1,MPI_integer,
     .    MPI_Max,MPI_Comm_World,MPIerror)
        ncmax = ntmp
#else
        nbasisloc = nbasis
        nbasismax = nbasis
#endif
        if (nbasisloc .gt. maxoloc)  overflow = .true.
        if (ncmax .gt. maxnc)  overflow = .true.

C  If overflow re-dimension arrays and repeat previous step
        if (overflow) then
          maxnc = ncmax
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
          allocate(listcold(maxnc,maxoloc))
          call memory('A','I',maxnc*maxoloc,'ordern')
          allocate(c(maxnc,maxoloc))
          call memory('A','D',maxnc*maxoloc,'ordern')
          allocate(cold(maxnc,maxoloc))
          call memory('A','D',maxnc*maxoloc,'ordern')
          overflow = .false.
          goto 10
        endif

        nbandsloc = ((nbands-1)/Nodes) + 1

        if (ioptlwf .ne. 0) then
   20     continue

C   This call was modified by DSP, Aug. 1998.
          call cspa(ioptlwf,iopt,natoms,nbasis,nbasisloc,lasto,isa,
     .        qa,enum,rcoor,rh,cell,xa,maxnh,numh,listhptr,listh,
     .        maxnc,c,numc,listc,ncmax,nctmax,nbands,overflow)

          nbandsloc = ((nbands-1)/Nodes) + 1

C         Check that dimensions are large enough 
          if (nbasis   .gt. maxo)    overflow = .true.
          if (nbasismax.gt. maxoloc) overflow = .true.
          if (nbandsloc.gt. maxlwf)  overflow = .true.
          if (ncmax    .gt. maxnc)   overflow = .true.
   
C Resize arrays while preserving contents of old version
          if (overflow) then
            if (nbasis .gt. maxo) then
              call memory('D','I',size(numc),'ordern')
              deallocate(numc)
              allocate(numc(nbasis))
              call memory('A','I',nbasis,'ordern')
C
              allocate(ntmp1(maxo))
              call memory('A','I',maxo,'ordern')
C
              ntmp1(1:maxo) = numcold(1:maxo)
              call memory('D','I',size(numcold),'ordern')
              deallocate(numcold)
              allocate(numcold(nbasis))
              call memory('A','I',nbasis,'ordern')
              numcold(1:maxo) = ntmp1(1:maxo)
C
              call memory('D','I',size(ntmp1),'ordern')
              deallocate(ntmp1)
            endif
            if (nbasis .gt. maxo .or. ncmax .gt. maxnc) then
              call memory('D','I',size(listc),'ordern')
              deallocate(listc)
              allocate(listc(ncmax,nbasis))
              call memory('A','I',ncmax*nbasis,'ordern')
            endif
            if (nbasismax .gt. maxoloc .or. ncmax .gt. maxnc) then
              allocate(rtmp2(maxnc,maxoloc))
              call memory('A','D',maxnc*maxoloc,'ordern')
C
              rtmp2(1:maxnc,1:maxoloc) = c(1:maxnc,1:maxoloc)
              call memory('D','D',size(c),'ordern')
              deallocate(c) 
              allocate(c(ncmax,nbasismax))
              call memory('A','D',ncmax*nbasismax,'ordern')
              c(1:maxnc,1:maxoloc) = rtmp2(1:maxnc,1:maxoloc)
C             
              rtmp2(1:maxnc,1:maxoloc) = cold(1:maxnc,1:maxoloc)
              call memory('D','D',size(cold),'ordern')
              deallocate(cold)
              allocate(cold(ncmax,nbasismax))
              call memory('A','D',ncmax*nbasismax,'ordern')
              cold(1:maxnc,1:maxoloc) = rtmp2(1:maxnc,1:maxoloc)
C             
              call memory('D','D',size(rtmp2),'ordern')
              deallocate(rtmp2)
              allocate(ntmp2(maxnc,maxoloc))
              call memory('A','I',maxnc*maxoloc,'ordern')
C             
              ntmp2(1:maxnc,1:maxoloc) = listcold(1:maxnc,1:maxoloc)
              call memory('D','I',size(listcold),'ordern')
              deallocate(listcold)
              allocate(listcold(ncmax,nbasis))
              call memory('A','I',ncmax*nbasis,'ordern')
              listcold(1:maxnc,1:maxoloc) = ntmp2(1:maxnc,1:maxoloc)
C             
              call memory('D','I',size(ntmp2),'ordern')
              deallocate(ntmp2)
            endif
            maxo = nbasis
            maxoloc = nbasismax
            maxlwf = nbandsloc
            maxnc = max(ncmax,maxnc)
            overflow = .false.
            goto 20
          endif
          maxnct = nctmax

        else
          
          nbandsloc = ((nbands-1)/Nodes) + 1

        endif
      endif
C ..................

C Calculate Chemical Potential, Max and Min eigenvalues, energy gap
C and HOMO and LUMO levels ..........

      if (chebef) then

        call timer( 'chempot', 1 )
        call chempot(h,s,numh,listhptr,listh,rcoorcp,ipcheb,beta,lasto,
     .               cell,xa,enum,nbasis,nbasisloc,natoms,maxnh,
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
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            call WhichNodeOrb(io,Nodes,BNode)
            if (Node.eq.BNode) then
              do in = 1,numcold(io)
                listc(in,io) = listcold(in,iio)
              enddo
            endif
#ifdef MPI
            call MPI_Bcast(listc(1,io),numcold(io),MPI_integer,
     .        BNode,MPI_Comm_World,MPIerror)
#endif
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
        call extrapolon(istep,iord,nspin,nbasisloc,maxo,maxoloc,
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
      call cgwf(iscf,itmax,ftol,eta,enum,nbasis,nbands,maxnh,
     .          numh,listhptr,listh,maxnc,numc,listc,h,s,c,
     .          g,hg,xi,fe,iter,dm,edm,nbasisloc,maxnct)
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
      nh = listhptr(nbasisloc) + numh(nbasisloc)
      qtot = 0.d0
      do io = 1,nh
        qtot = qtot + dm(io) * s(io)
      enddo

#ifdef MPI
C Global reduction of qtot
#ifdef NODAT
      call MPI_AllReduce(qtot,rtmp,1,MPI_double_precision,MPI_sum,
     .  MPI_Comm_World,MPIerror)
#else
      call MPI_AllReduce(qtot,rtmp,1,DAT_double,MPI_sum,
     .  MPI_Comm_World,MPIerror)
#endif
      qtot = rtmp
#endif

      Ecorrec = eta * (enum - qtot)
      if (Node.eq.0) then
      write(6,"(a,f12.4)") 'ordern: qtot (after  DM normalization) = ',
     .   qtot
      endif
C ..................

C Save LWF's to disk .................................................
      do io = 1,nbasisloc
        call LocalToGlobalOrb(io,Node,Nodes,iio)
        do in = 1,numcold(iio)
          listcold(in,io) = listc(in,iio)
        enddo
      enddo
      if (Node.eq.0) then
        writelwf = fdf_boolean('WriteLWF', .true. )
      endif
#ifdef MPI
      call MPI_Bcast(writelwf,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if (writelwf) then
        call iolwf( 'write', maxnc, ncmax, maxoloc, nbasis, 1,
     .             numc, listcold, c, found )
      endif
C ....................

*     call timer( 'ordern', 2 )

      end subroutine ordern

      end module m_ordern
