module m_pexsi
  public :: pexsi
CONTAINS
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
    integer :: MPIerror, stat(MPI_STATUS_SIZE)
    integer :: block_size, norbs_slack, nnz_slack
    integer :: ispin, maxnhtot, ih

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
call MPI_AllReduce(maxnh,maxnhtot,1,MPI_integer,MPI_sum,MPI_Comm_World,MPIerror)
nnz   = maxnhtot

! Compute the number of orbs that we must hand over to the
! last processor from the first
block_size = no_u/mpisize
norbs_slack = no_u - mpisize*block_size

if (mpirank == 0) then

  numColLocal = block_size 
  nnzLocal = sum(numh(1:block_size))  ! discard the rest
  nnz_slack = sum(numh(block_size+1:no_l))  
  ! send the information to the last processor
  call MPI_Send(numh(block_size+1),norbs_slack,MPI_Integer,mpisize-1,0,MPI_comm,ierr)
  print *, "First node: no_l, numColLocal, maxnh, nnzLocal: ", no_l, numColLocal, maxnh, nnzLocal
  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo

else if (mpirank == mpisize-1) then

  numColLocal = block_size + norbs_slack
  allocate(tmpi(norbs_slack))
  call MPI_Recv(tmpi,norbs_slack,MPI_Integer,0,0,MPI_comm,stat,ierr)
  nnz_slack = sum(tmpi(1:norbs_slack))  
  nnzLocal = sum(numh(1:block_size)) + nnz_slack
  print *, "Last node: no_l, numColLocal, maxnh, nnzLocal: ", no_l, numColLocal, maxnh, nnzLocal

  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,block_size
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo
  do ih = block_size+1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + tmpi(ih)
  enddo

else

  numColLocal = block_size 
  nnzLocal = sum(numh(1:block_size))
  print *, "Node ", mpirank, ": no_l, numColLocal, maxnh, nnzLocal: ", no_l, numColLocal, maxnh, nnzLocal

  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo

endif

! Prepare the matrix indexes and values
if (mpirank == 0) then
  call MPI_Send(listh(nnzLocal+1),nnz_slack,MPI_Integer,mpisize-1,1,MPI_comm,ierr)
  call MPI_Send(S(nnzLocal+1),nnz_slack,MPI_Double_Precision,mpisize-1,2,MPI_comm,ierr)
  call MPI_Send(H(nnzLocal+1,ispin),nnz_slack,MPI_Double_Precision,mpisize-1,3,MPI_comm,ierr)
   rowindLocal => listh
   HnzvalLocal => H(:,ispin)
   SnzvalLocal => S
   DMnzvalLocal => DM(:,ispin)
   EDMnzvalLocal => EDM(:,ispin)
   allocate(FDMnzvalLocal(1:maxnh))  ! Note expanded array for future use

else if (mpirank == mpisize-1) then

   allocate(rowindLocal(1:nnzLocal))
   allocate(HnzvalLocal(1:nnzLocal))
   allocate(SnzvalLocal(1:nnzLocal))
   allocate(DMnzvalLocal(1:nnzLocal))
   allocate(EDMnzvalLocal(1:nnzLocal))
   allocate(FDMnzvalLocal(1:nnzLocal))

   rowindLocal(1:maxnh) =  listh(1:maxnh)
   HnzvalLocal(1:maxnh) =  H(1:maxnh,ispin)
   SnzvalLocal(1:maxnh) =  S(1:maxnh)
   call MPI_Recv(rowindLocal(maxnh+1),nnz_slack,MPI_Integer,0,1,MPI_comm,stat,ierr)
   call MPI_Recv(SnzvalLocal(maxnh+1),nnz_slack,MPI_Double_Precision,0,2,MPI_comm,stat,ierr)
   call MPI_Recv(HnzvalLocal(maxnh+1),nnz_slack,MPI_Double_Precision,0,3,MPI_comm,stat,ierr)

else
   rowindLocal => listh
   HnzvalLocal => H(:,ispin)
   SnzvalLocal => S
   DMnzvalLocal => DM(:,ispin)
   EDMnzvalLocal => EDM(:,ispin)
   allocate(FDMnzvalLocal(1:nnzLocal))

endif

call MPI_Barrier(MPI_comm,ierr)

   
! Data is for the DNA matrix.
temperature      = 300.0d0    ! Units??
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


  DMnzvalLocal(:) = 2*HnzvalLocal(:)
  EDMnzvalLocal(:) = -1.0_dp

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


call MPI_Barrier(MPI_comm,ierr)

! Recover DM and EDM
if (mpirank == 0) then

 DM(1:nnzLocal,ispin) = DMnzvalLocal(1:nnzLocal)
 EDM(1:nnzLocal,ispin) = EDMnzvalLocal(1:nnzLocal)
 call MPI_Recv(DM(nnzLocal+1,ispin),nnz_slack,MPI_Double_Precision,mpisize-1,1,MPI_comm,stat,ierr)
 call MPI_Recv(EDM(nnzLocal+1,ispin),nnz_slack,MPI_Double_Precision,mpisize-1,2,MPI_comm,stat,ierr)

 print *, "Node: ", mpirank, ": H: ",  H(nnzLocal+1:nnzLocal+2,ispin)
 print *, "Node: ", mpirank, ": DM: ", DM(nnzLocal+1:nnzLocal+2,ispin)
 print *, "Node: ", mpirank, ": EDM: ", EDM(nnzLocal+1:nnzLocal+2,ispin)

else if (mpirank == mpisize-1) then

   call MPI_Send(DMnzvalLocal(maxnh+1),nnz_slack,MPI_Double_Precision,0,1,MPI_comm,ierr)
   call MPI_Send(EDMnzvalLocal(maxnh+1),nnz_slack,MPI_Double_Precision,0,2,MPI_comm,ierr)
   DM(1:maxnh,ispin) = DMnzvalLocal(1:maxnh)
   EDM(1:maxnh,ispin) = EDMnzvalLocal(1:maxnh)

   deallocate(rowindLocal)
   deallocate(HnzvalLocal)
   deallocate(SnzvalLocal)
   deallocate(DMnzvalLocal)
   deallocate(EDMnzvalLocal)
   deallocate(FDMnzvalLocal)

else

   deallocate(FDMnzvalLocal)

endif

deallocate(colPtrLocal)

    call timer("pexsi", 2)
#endif 

  end subroutine pexsi
end module m_pexsi
