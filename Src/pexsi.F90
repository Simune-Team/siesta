module m_pexsi

  public :: pexsi

CONTAINS

! SIESTA interface to PEXSI
!
! This version uses the standard PEXSI distribution (by tricking Siesta into
! using it)
!
  subroutine pexsi(no_u, no_l, nspin,  &
       maxnh, numh, listhptr, listh, H, S, qtot, DM, EDM)

    use precision, only  : dp
#ifdef MPI
    use mpi_siesta
#endif
    implicit          none

    integer, intent(in)  :: maxnh, no_u, no_l, nspin
    integer, intent(in), target  :: listh(maxnh), numh(*), listhptr(*)
    real(dp), intent(in), target :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: qtot
    real(dp), intent(out), target:: DM(maxnh,nspin), EDM(maxnh,nspin)

#ifndef MPI
    call die("PEXSI needs MPI")
#else

    integer :: mpi_comm = mpi_comm_world
    integer :: MPIerror, stat(MPI_STATUS_SIZE), count
    integer :: bs, norbs_slack, nnz_slack
    integer :: ispin, maxnhtot, ih, nnzold

!Lin variables
integer :: nrows, nnz, nnzLocal, numColLocal
integer, pointer, dimension(:) ::  tmpi => null()
integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
!
real(dp), pointer, dimension(:) :: &
	HnzvalLocal=>null(), SnzvalLocal=>null(),  &
        DMnzvalLocal => null() , EDMnzvalLocal => null(), &
	FDMnzvalLocal => null()
!
integer :: numPole
real(dp) :: temperature, numElectronExact, numElectron,&
	gap, deltaE
real(dp) :: mu, muMin, muMax
integer:: muMaxIter
real(dp) :: poleTolerance, numElectronTolerance
integer:: npPerPole, nprow, npcol
integer :: mpirank, mpisize, ierr
!------------


external      :: timer

call timer("pexsi", 1)

ispin = 1
if (nspin /=1) then
   call die("Spin polarization not yet supported in PEXSI")
endif

call mpi_comm_rank( MPI_COMM, mpirank, ierr )
call mpi_comm_size( MPI_COMM, mpisize, ierr )

call MPI_Barrier(MPI_comm,ierr)

nrows = no_u
numColLocal = no_l

nnzLocal = sum(numh(1:numColLocal))
call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,MPI_Comm,MPIerror)


!  print *, "Node ", mpirank, ": no_l, numColLocal, maxnh, nnzLocal: ", &
!           no_l, numColLocal, maxnh, nnzLocal

  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo

   allocate(rowindLocal(1:nnzLocal))
   allocate(HnzvalLocal(1:nnzLocal))
   allocate(SnzvalLocal(1:nnzLocal))
   allocate(DMnzvalLocal(1:nnzLocal))
   allocate(EDMnzvalLocal(1:nnzLocal))
   allocate(FDMnzvalLocal(1:nnzLocal))

   rowindLocal(1:nnzLocal) = listh(1:nnzLocal)
   HnzvalLocal(1:nnzLocal) = H(1:nnzLocal,ispin)
   SnzvalLocal(1:nnzLocal) = S(1:nnzLocal)

! Data is for the DNA matrix.
temperature      = 3000.0d0    ! Units??
numElectronExact = qtot   ! 2442.0d0 for DNA
numPole          = 20
gap              = 0.0d0
! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
deltaE           = 3.0d0
! Initial guess of chemical potential, also updated after pexsi.
mu               = -0.60d0
! Lower/Upper bound for the chemical potential.
muMin            = -1.0d0
muMax            =  0.0d0
! muMaxIter should be 1 or 2 later when combined with SCF.
muMaxIter        = 10
! Do not compute a pole if the corresponding weight is < poleTolerance.
poleTolerance    = 1d-8
! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
numElectronTolerance = 1d-1
! Number of processors used for each pole. At the moment use mpisize.
! Later can be changed to 
npPerPole        = mpisize


!!$call f_ppexsi_interface( &
!!$	nrows,&
!!$	nnz,&
!!$	nnzLocal,&
!!$	numColLocal,&
!!$	colptrLocal,&
!!$        rowindLocal,&
!!$	HnzvalLocal,&
!!$	SnzvalLocal,&
!!$	DMnzvalLocal,&
!!$	EDMnzvalLocal,&
!!$	FDMnzvalLocal,&
!!$	numPole,&
!!$	temperature,&
!!$	numElectronExact,&
!!$	numElectron,&
!!$	gap,&
!!$	deltaE,&
!!$	mu,&
!!$	muMin,&
!!$        muMax,&
!!$	muMaxIter,&
!!$	poleTolerance,&
!!$	numElectronTolerance,&
!!$	MPI_COMM,&
!!$	npPerPole )

if( mpirank == 0 ) then
	write(*, *) "mu          = ", mu
	write(*, *) "numElectron = ", numElectron
endif

   DM(1:nnzLocal,ispin) = DMnzvalLocal(1:nnzLocal)
   EDM(1:nnzLocal,ispin) = EDMnzvalLocal(1:nnzLocal)

  deallocate(rowindLocal)
  deallocate(HnzvalLocal)
  deallocate(SnzvalLocal)
  deallocate(DMnzvalLocal)
  deallocate(EDMnzvalLocal)
  deallocate(FDMnzvalLocal)
  deallocate(colPtrLocal)

    call timer("pexsi", 2)
#endif 

  end subroutine pexsi
end module m_pexsi
