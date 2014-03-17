      module m_evolve
      
      CONTAINS
      subroutine evolve(no, nspin, maxspn, maxuo, maxnh, maxnd,
     .                  maxo, numh, listhptr, listh, numd,
     .                  listdptr, listd, H, S, eo,
     .                  gamma, xij, indxuo, nk, kpoint, wk,
     .                  Dnew, Enew, nuotot, delt)

C *********************************************************************
C Subroutine to time-evolve the eigenvectors, calculate the density 
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization).
C Written by A. Tsolakidis, May 2000 after a suboutine
C by P. Ordejon and J. M. Soler.
C Modified by D. Sanchez-Portal, November 2002
C Modified by D. Sanchez-Portal, March 2008
C Modified by M. Ahsan Zeb, 11 March 2011.
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals the supercell
C integer nspin               : Spin polarization (1 or 2)
C integer maxspn              : Maximum number of spin orentations 
C integer maxuo               : Maximum number of orbitals stored in a 
C                               given Node
C integer maxnh               : Maximum number of orbitals interacting
C integer maxnd               : Maximum number of nonzero elements of
C                               each row of density matrix
C integer maxo                : Maximum number of orbitals in the unit cell
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix
C integer listh(maxlh)        : Nonzero hamiltonian-matrix element
C                               column indexes for each matrix row
C integer numd(nuo)           : Number of nonzero elements of each row
C                               of density matrix
C integer listdptr(nuo)       : Pointer to each row (-1) of the
C                               density matrix
C integer listd(maxnh)        : Nonzero density-matrix element column
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C logical gamma               : Only gamma point?
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer nk                  : Number of k points
C real*8  kpoint(3,nk)        : k point vectors
C real*8  wk(nk)              : k point weights (must sum one)
C integer nuotot              : total number of orbitals in unit cell
C                               over all processors
C double precision delt       : time step in units of the inverse of 
C                               the Hamiltonian
C *************************** OUTPUT **********************************
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix
C real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
C real*8 eo(maxo,maxspn,nk)   : Output instantaneous eigenvalues 
C                              (only calculated if explicitly required 
C                               by user and in the last MD step) 
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C Enew returned in the units of H.
C delt in femtoseconds
C *************************** Parallel ********************************
C  Very important!!!!!!!!!!!!
C  DSP: This subroutine is not yet prepared to run in parallel!!!!
C  Sorry!
C When running in parallel some of the dimensions are now the
C maximum per node and the corresponding number passed in as
C an argument is the number of locally stored values. The
C variables for which this is the case are:
C
C maxuo/no
C
C *********************************************************************
C
C  Modules
C
      use precision
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs
      use fdf
      use densematrix,  only : Haux, Saux
      use alloc
      use m_memory
      use sys, only: die
#ifdef MPI
      use mpi_siesta,   only : mpi_bcast, mpi_comm_world,
     $                         mpi_logical
#endif


      implicit none

      integer
     .  maxnd, maxnh, maxspn, maxuo, maxo, nk, no, nspin, nuotot

      integer
     .  indxuo(no), listh(maxnh), numh(*), listd(maxnd), numd(*),
     .  listhptr(*), listdptr(*)

      double precision
     .  Dnew(maxnd,nspin), Enew(maxnd,nspin),
     .  H(maxnh,nspin), kpoint(3,nk),
     .  S(maxnh), wk(nk),
     .  xij(3,maxnh), delt, eo(maxo,maxspn,nk)

      logical
     .  gamma

      external
     . io_assign, io_close

#ifdef MPI
      integer
     .  MPIerror
      external
     .  diagkp
#endif
      logical
     .  frstme

C  Internal variables .............................................
      integer
     .  io, iuo, iu, naux, nhs,  nuo

C Dynamic arrays
      real(dp), dimension(:), allocatable, save :: aux,aux2
      real(dp), pointer, save :: Dk(:)
      real(dp), pointer, save :: Ek(:)

#ifdef MPI
      logical, save ::
     .  ParallelOverK
#endif
C  ....................

      save frstme
      data frstme /.true./


C Get Node number and calculate local orbital range
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
      if (frstme) then
        if (Node.eq.0) then
          ParallelOverK = fdf_boolean( 'Diag.ParallelOverK', .false. )
        endif
        call MPI_Bcast(ParallelOverK,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
      endif
#else
      Node = 0
      Nodes = 1
      nuo = nuotot
#endif
C  ....................

C Start time counter ................................................
      call timer( 'evolve', 1 )
C ......................

C Check internal dimensions ..........................................
        nhs  =  2 * nuotot * nuo
        naux = 2 * nuotot        
#ifdef MPI
        if (ParallelOverK) then
          nhs  = 2 * nuotot * nuotot
        endif
#endif

C Allocate local arrays
      call re_alloc(Haux,1,nhs,name='Haux',routine='evolve')
      call re_alloc(Saux,1,nhs,name='Saux',routine='evolve')
      if(frstme) then 
        if(.not.gamma) then 
          allocate(Dk(nhs),stat=mem_stat)
          call memory('A','D',nhs,'evolve',stat=mem_stat)
          allocate(Ek(nhs),stat=mem_stat)
          call memory('A','D',nhs,'evolve',stat=mem_stat)
        endif
        allocate(aux(naux),stat=mem_stat)
        call memory('A','D',naux,'evolve',stat=mem_stat)
        allocate(aux2(naux),stat=mem_stat)
        call memory('A','D',naux,'evolve',stat=mem_stat)
        frstme=.false.
      endif
C.....................................................................

C DSP
C Here diagon checks for the consistency of the orbital 
C indices in different nodes. Since this is not yet
C desinged to run in parallel, I skip this check

C Call apropriate routine .............................................

      if (nspin.le.2 .and. gamma) then
        call evolg( nspin, nuo, no, maxo, maxnh, maxnd,
     .              numh, listhptr, listh, numd,
     .              listdptr, listd, H, S, eo,
     .              Dnew, Enew, Haux, Saux,
     .              nuotot, delt)
      elseif (nspin.le.2 .and. .not.gamma) then
        call evolk( nspin, maxspn, nuo, no, maxo, maxnh, maxnd,
     .              numh, listhptr, listh, numd, listdptr,
     .              listd, H, S, eo,
     .              xij, indxuo, nk, kpoint, wk,
     .              Dnew, Enew, 
     .              Haux, Saux, Dk, Ek, 
     .              nuotot, delt)

      elseif (nspin.eq.4 .and. gamma) then 
         stop 'evolve: ERROR: non-collinear spin not yet implemented'
      elseif (nspin.eq.4 .and. .not. gamma) then 
         stop 'evolve: ERROR: non-collinear spin not yet implemented'
      endif 
C ....................


C Stop time counter ...................................................
      call timer( 'evolve', 2 )
C .......................

      end subroutine evolve
      
      END module m_evolve










