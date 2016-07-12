      module m_evolve
      
      CONTAINS
      subroutine evolve(no, nspin, maxspn, maxuo, maxnh, maxnd,        &
                       maxo, gamma, indxuo, nk, kpoint, wk,            &
                       Dnew, Enew, nuotot, delt, istep, itded)         

! *********************************************************************
! Subroutine to time-evolve the eigenvectors, calculate the density 
! and energy-density matrices, and occupation weights of each 
! eigenvector, for given Hamiltonian and Overlap matrices (including
! spin polarization).
! Written by A. Tsolakidis, May 2000 after a suboutine
! by P. Ordejon and J. M. Soler.
! Modified by D. Sanchez-Portal, November 2002
! Modified by D. Sanchez-Portal, March 2008
! Modified by Rafi Ullah, October ,2015 making it parallel using 
! Matrix Swtich.
! **************************** INPUT **********************************
! integer no                  : Number of basis orbitals the supercell
! integer nspin               : Spin polarization (1 or 2)
! integer maxspn              : Maximum number of spin orentations 
! integer maxuo               : Maximum number of orbitals stored in a 
!                               given Node
! integer maxnh               : Maximum number of orbitals interacting
! integer maxnd               : Maximum number of nonzero elements of
!                               each row of density matrix
! integer maxo                : Maximum number of orbitals in the unit cell
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to each row (-1) of the
!                               hamiltonian matrix
! integer listh(maxlh)        : Nonzero hamiltonian-matrix element
!                               column indexes for each matrix row
! integer numd(nuo)           : Number of nonzero elements of each row
!                               of density matrix
! integer listdptr(nuo)       : Pointer to each row (-1) of the
!                               density matrix
! integer listd(maxnh)        : Nonzero density-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! logical gamma               : Only gamma point?
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not used if only gamma point)
! integer indxuo(no)          : Index of equivalent orbital in unit cell
!                               Unit cell orbitals must be the first in
!                               orbital lists, i.e. indxuo.le.nuo, with
!                               nuo the number of orbitals in unit cell
! integer nk                  : Number of k points
! real*8  kpoint(3,nk)        : k point vectors
! real*8  wk(nk)              : k point weights (must sum one)
! integer nuotot              : total number of orbitals in unit cell
!                               over all processors
! double precision delt       : time step in units of the inverse of 
!                               the Hamiltonian
! *************************** OUTPUT **********************************
! real*8 Dnew(maxnd,nspin)    : Output Density Matrix
! real*8 Enew(maxnd,nspin)    : Output Energy-Density Matrix
! real*8 eo(maxo,maxspn,nk)   : Output instantaneous eigenvalues 
!                              (only calculated if explicitly required 
!                               by user and in the last MD step) 
! *************************** UNITS ***********************************
! xij and kpoint must be in reciprocal coordinates of each other.
! Enew returned in the units of H.
! delt in femtoseconds
! *************************** Parallel ********************************
! When running in parallel some of the dimensions are now the
! maximum per node and the corresponding number passed in as
! an argument is the number of locally stored values. The
! variables for which this is the case are:
! 
! maxuo/no
! 
! *********************************************************************
! 
!  Modules
  
      use precision
      use parallel,          only : Node, Nodes, BlockSize
      use parallelsubs,      only : GlobalToLocalOrb, GetNodeOrbs
      use fdf
      use alloc
      use m_memory
!      use densematrix,      only : Haux, Saux,psi
!      use sparse_matrices,  only : S, H, numh, listh, listhptr, xijo 
!      use m_eo,             only : eo
      use sys,              only : die
      use MatrixSwitch
#ifdef MPI
      use mpi_siesta,       only : MPI_Bcast, MPI_Comm_World,MPI_logical
#endif
      !
      implicit none
      !
      integer, intent(in)      ::  maxnd, maxnh, maxspn, maxuo, maxo, nk, no, nspin 
      integer                  ::  nuotot, istep, itded,  indxuo(no) 
      double precision         ::  Dnew(maxnd,nspin), Enew(maxnd,nspin) 
      double precision         ::  kpoint(3,nk), wk(nk), delt
      logical                  ::  gamma
      external                 ::  io_assign, io_close
      !
#ifdef MPI
      integer                  ::  MPIerror
      external                 ::  diagkp
#endif
      logical, save            ::  frstme  = .true.
      ! Internal variables ...
      integer                  ::  io, iuo, iu, naux, nhs,  nuo,npsi
      real(dp), pointer, save  :: Dk(:)
      real(dp), pointer, save  :: Ek(:)
      real(dp), dimension(:), allocatable, save :: aux,aux2
#ifdef MPI
      logical, save            :: ParallelOverK
#endif
     !
#ifdef MPI
      call GetNodeOrbs(nuotot,Node,Nodes,nuo)
      if (frstme) then
        if (Node.eq.0) then
          ParallelOverK = fdf_boolean( 'Diag.ParallelOverK', .false. )
        endif
        if (ParallelOverK) then
          call MPI_Bcast(ParallelOverK,1,MPI_logical,0,MPI_Comm_World, MPIerror)
        end if
      endif
#else
      Node = 0
      Nodes = 1
      nuo = nuotot
#endif
      ! Start time counter ................................................
      call timer( 'evolve', 1 )
      nhs  =  2 * nuotot * 2* nuotot
      naux = 2 * nuotot        
      npsi=nuotot*maxuo*nspin
#ifdef MPI
      if (ParallelOverK) then
        nhs  = 2 * nuotot * nuotot
      endif
#endif
      !
!      call re_alloc(psi,1,npsi,name='psi',routine='evolve')
!      call re_alloc(psi2,1,npsi,name='psi2',routine='evolve')
!      call re_alloc(Haux,1,nhs,name='Haux',routine='evolve')
!      call re_alloc(Saux,1,nhs,name='Saux',routine='evolve')
      !
!      if(frstme) then 
!        if(.not.gamma) then 
!          allocate(Dk(nhs),stat=mem_stat)
!          call memory('A','D',nhs,'evolve',stat=mem_stat)
!          allocate(Ek(nhs),stat=mem_stat)
!          call memory('A','D',nhs,'evolve',stat=mem_stat)
!        endif
!        allocate(aux(naux),stat=mem_stat)
!        call memory('A','D',naux,'evolve',stat=mem_stat)
!        allocate(aux2(naux),stat=mem_stat)
!        call memory('A','D',naux,'evolve',stat=mem_stat)
!        frstme=.false.
!      endif
      ! Call apropriate routine .............................................
      if (nspin.le.2 .and. gamma) then
        call evolg( nspin, nuo, no, maxo, maxnh, maxnd,                   &
                    Dnew, Enew, nuotot, delt,istep,itded)
      elseif (nspin.le.2 .and. .not.gamma) then
      stop 'evolve: Error: kpoint for TDED is not implimented yet'
 !       call evolk( nspin, maxspn, nuo, no, maxo, maxnh, maxnd,           &
 !                   indxuo, nk, kpoint, wk, Dnew, Enew, Dk, Ek,           &
 !                   nuotot, delt,Haux,Saux,psi)

      elseif (nspin.eq.4 .and. gamma) then 
        stop 'evolve: ERROR: non-collinear spin not yet implemented'
      elseif (nspin.eq.4 .and. .not. gamma) then 
        stop 'evolve: ERROR: non-collinear spin not yet implemented'
      endif 
      ! Stop time counter ...................................................
      call timer( 'evolve', 2 )
      !
      end subroutine evolve
      END module m_evolve
