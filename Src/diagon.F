      subroutine diagon(no, nspin, maxspn, maxuo, maxnh, maxnd, 
     .                  maxo, numh, listhptr, listh, numd, 
     .                  listdptr, listd, H, S,
     .                  qtot, fixspin, qs, temp, e1, e2,
     .                  gamma, xij, indxuo, nk, kpoint, wk,
     .                  eo, qo, Dnew, Enew, ef, efs, Entropy, nuotot)
C *********************************************************************
C Subroutine to calculate the eigenvalues and eigenvectors, density
C and energy-density matrices, and occupation weights of each 
C eigenvector, for given Hamiltonian and Overlap matrices (including
C spin polarization).
C Written by J.Soler and P.Ordejon, August 1997.
C Now includes the option to parallelise run over K points.
C Modified by J.Gale, November 1999.
C **************************** INPUT **********************************
C integer no                  : Number of basis orbitals
C integer nspin               : Spin polarization (1 or 2)
C integer maxspn              : Second dimension of eo and qo
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : Maximum number of nonzero elements of 
C                               each row of density matrix
C integer maxo                : First dimension of eo and qo
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
C real*8  qtot                : Number of electrons in unit cell
C logical fixspin             : Fix the spin of the system?
C real*8  qs(nspin)           : Number of electrons in unit cell for each
C                               spin component (if fixed spin option is used)
C real*8  temp                : Electronic temperature 
C real*8  e1, e2              : Energy range for density-matrix states
C                               (to find local density of states)
C                               Not used if e1 > e2
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
C *************************** OUTPUT **********************************
C real*8 eo(maxo,maxspn,nk)   : Eigenvalues
C real*8 qo(maxo,maxspn,nk)   : Occupations of eigenstates
C real*8 Dnew(maxnd,nspin)    : Output Density Matrix
C real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
C real*8 ef                   : Fermi energy
C real*8 efs(nspin)           : Fermi energy for each spin
C                                (for fixed spin calculations)
C real*8 Entropy              : Electronic entropy
C *************************** UNITS ***********************************
C xij and kpoint must be in reciprocal coordinates of each other.
C temp and H must be in the same energy units.
C eo, Enew and ef returned in the units of H.
C *************************** Parallel ********************************
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
      use parallel
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer
     .  maxnd, maxnh, maxspn, maxuo, maxo, nk, no, nspin, nuotot

      integer 
     .  indxuo(no), listh(maxnh), numh(*), listd(maxnd), numd(*),
     .  listhptr(*), listdptr(*)

      double precision
     .  Dnew(maxnd,nspin), e1, e2, ef, Enew(maxnd,nspin), Entropy,
     .  eo(maxo,maxspn,nk), H(maxnh,nspin), kpoint(3,nk), 
     .  qo(maxo,maxspn,nk), qtot, S(maxnh), temp, wk(nk), 
     .  xij(3,maxnh), qs(nspin), efs(nspin), qspiral(3)
     
      logical
     .  fixspin, gamma, spiral

      external
     .  diagg, diagk, diag2g, diag2k, io_assign, io_close, memory,
     .  readsp, diagsprl

#ifdef MPI
      integer
     .  MPIerror
      external
     .  diagkp
      logical
     .  frstme
#endif

      logical
     .  getD, getPSI

      integer
     .  io, iuo, iu, naux, nhs, npsi, nuo, Node, Nodes

C     Dynamic arrays
      double precision, dimension(:), allocatable, save :: Haux, Saux
      double precision, dimension(:), allocatable, save :: psi, aux
      integer         , dimension(:), allocatable, save :: muo

#ifdef MPI
      logical, save ::
     .  ParallelOverK
      save frstme
      data frstme /.true./
#endif
C  ....................

C Get Node number and calculate local orbital range
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
      if (frstme) then
        if (Node.eq.0) then
          ParallelOverK = fdf_boolean( 'ParallelOverK', .false. )
        endif
        call MPI_Bcast(ParallelOverK,1,MPI_logical,0,MPI_Comm_World,
     .    MPIerror)
        frstme = .false.
      endif
#else
      Node = 0
      Nodes = 1
      nuo = nuotot
#endif

C Start time counter ................................................
      call timer( 'diagon', 1 )
C ......................

C Read spin-spiral wavevector (if defined)..........................
      call readsp( qspiral, spiral )
      if (spiral) then
        if (gamma) write(6,*) 
     .    'diagon: WARNING: spin-spiral requires k sampling'
        if (nspin.ne.4) write(6,*)
     .    'diagon: WARNING: spin-spiral requires nspin=4'
      endif
C ......................

C Check internal dimensions ..........................................
      if (nspin.le.2 .and. gamma) then
        nhs  = nuotot * nuo
        npsi = nuotot * maxuo * nspin
        naux = nuotot * 5
      elseif (nspin.le.2 .and. .not.gamma) then
        nhs  = 2 * nuotot * nuo
        npsi = 2 * nuotot * nuo
        naux = 2 * nuotot * 5            !!!! AG: Put factor 2 !!
#ifdef MPI
        if (ParallelOverK) then
          nhs  = 2 * nuotot * nuotot
          npsi = 2 * nuotot * nuotot
        endif
#endif
      elseif (nspin.eq.4) then
        nhs  = 2 * (2*nuotot) * (2*nuo)
        npsi = 2 * (2*nuotot) * (2*maxuo)
        naux = (2*nuotot) * 5
      else
        if (Node.eq.0) then
          stop 'diagon: ERROR: incorrect value of nspin'
        else
          stop
        endif
      endif

C Allocate local arrays
      allocate(Haux(nhs))
      call memory('A','D',nhs,'diagon')
      allocate(Saux(nhs))
      call memory('A','D',nhs,'diagon')
      allocate(psi(npsi))
      call memory('A','D',npsi,'diagon')
      allocate(aux(naux))
      call memory('A','D',naux,'diagon')
      allocate(muo(nuo))
      call memory('A','I',nuo,'diagon')

C Check indxuo .......................................................
      do iuo = 1,nuo
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (iuo.le.0 .or. iuo.gt.nuotot) then
          if (Node.eq.0) then
            write(6,*) 'diagon: ERROR: invalid index: io, indxuo =',
     .        io, indxuo(io)
            stop 'diagon: ERROR: invalid indxuo'
          else
            stop
          endif
        endif
        call GlobalToLocalOrb(indxuo(io),Node,Nodes,iuo)
        if (iuo.gt.0) then
          muo(iuo) = muo(iuo) + 1
        endif
      enddo
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
          if (Node.eq.0) then
            write(6,'(/,2a,3i6)') 'diagon: ERROR: inconsistent indxuo',
     .       '. iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
            stop 'diagon: ERROR: inconsistent indxuo.'
          else
            stop
          endif
        endif
      enddo
C ......................

C Call apropriate routine .............................................
      getD = .true.
      getPSI = .false.

      if (nspin.le.2 .and. gamma) then
        call diagg( nspin, nuo, no, maxuo, maxnh, maxnd, maxo,
     .              numh, listhptr, listh, numd, listdptr, listd, 
     .              H, S, getD, getPSI, fixspin, qtot, qs, temp, e1, e2,
     .              eo, qo, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, aux, 
     .              nuotot, Node, Nodes)
      elseif (nspin.le.2 .and. .not.gamma) then
#ifdef MPI
        if (ParallelOverK) then
          call diagkp( nspin, nuo, no, maxspn, maxuo, maxnh, maxnd, 
     .              maxo, numh, listhptr, listh, numd, listdptr, 
     .              listd, H, S, getD, fixspin, qtot, qs, temp, 
     .              e1, e2, xij, indxuo, nk, kpoint, wk,
     .              eo, qo, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, Haux, Saux, aux, 
     .              nuotot, Node, Nodes)
        else
#endif
          call diagk( nspin, nuo, no, maxspn, maxuo, maxnh, maxnd, 
     .              maxo, numh, listhptr, listh, numd, listdptr,
     .              listd, H, S, getD, getPSI, fixspin, qtot, qs, temp, 
     .              e1, e2, xij, indxuo, nk, kpoint, wk,
     .              eo, qo, Dnew, Enew, ef, efs, Entropy,
     .              Haux, Saux, psi, Haux, Saux, aux, 
     .              nuotot, Node, Nodes)
#ifdef MPI
        endif
#endif
      elseif (nspin.eq.4 .and. gamma) then
        call diag2g( nuo, no, maxuo, maxnh, maxnd, maxo,
     .               numh, listhptr, listh, numd, listdptr,
     .               listd, H, S, getD, qtot, temp, e1, e2,
     .               eo, qo, Dnew, Enew, ef, Entropy,
     .               Haux, Saux, psi, aux, 
     .               nuotot, Node, Nodes)
      elseif (nspin.eq.4 .and. spiral) then
        call diagsprl( nuo, no, maxuo, maxnh, maxnd, maxo,
     .               numh, listhptr, listh, numd, listdptr,
     .               listd, H, S, getD, qtot, temp, e1, e2,
     .               xij, indxuo, nk, kpoint, wk,
     .               eo, qo, Dnew, Enew, ef, Entropy, qspiral,
     .               Haux, Saux, psi, Haux, Saux, aux,
     .               nuotot, Node, Nodes)
      elseif (nspin.eq.4 .and. .not.gamma) then
        call diag2k( nuo, no, maxuo, maxnh, maxnd, maxo,
     .               numh, listhptr, listh, numd, listdptr,
     .               listd, H, S, getD, qtot, temp, e1, e2,
     .               xij, indxuo, nk, kpoint, wk,
     .               eo, qo, Dnew, Enew, ef, Entropy,
     .               Haux, Saux, psi, Haux, Saux, aux, 
     .               nuotot, Node, Nodes)
      endif
C ....................

C Free local arrays
      call memory('D','I',size(muo),'diagon')
      deallocate(muo)
      call memory('D','D',size(aux),'diagon')
      deallocate(aux)
      call memory('D','D',size(psi),'diagon')
      deallocate(psi)
      call memory('D','D',size(Saux),'diagon')
      deallocate(Saux)
      call memory('D','D',size(Haux),'diagon')
      deallocate(Haux)

C Stop time counter ...................................................
      call timer( 'diagon', 2 )
C .......................

      end
