module m_pexsi_exchange

  public :: pexsi_exchange

CONTAINS

! SIESTA interface to PEXSI
!
! This version uses the standard block-cyclic SIESTA distribution, but with
! a blocksize of bs=Norbs/Nodes. Note: remainder=Norbs-bs*Nodes.
! This results in the first node holding bs+remainder orbitals, and the rest
! holding bs orbitals.
! In PEXSI, the last orbital is supposed to hold the slack (those 'remainder'
! last orbitals). Thus it is necessary to exchange information between the first
! and last orbitals.
!
  subroutine pexsi_EXCHANGE(no_u, no_l, nspin,  &
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
call MPI_AllReduce(maxnh,maxnhtot,1,MPI_integer,MPI_sum,MPI_Comm,MPIerror)
nnz   = maxnhtot

! Compute the number of orbs that we must hand over to the
! last processor from the first
bs = no_u/mpisize
norbs_slack = no_u - mpisize*bs

if (mpirank == 0) then

  numColLocal = bs 
  nnzLocal = sum(numh(1:bs))  ! discard the rest
  if (norbs_slack /= (no_l - bs)) call die("norbs_slack /= no_l-bs in 0")

  nnz_slack = sum(numh(bs+1:no_l))  
  ! send the information to the last processor
  call MPI_Send(numh(bs+1),norbs_slack,MPI_Integer, &
                mpisize-1,0,MPI_comm,ierr)
  print *, "First node: no_l, numColLocal, maxnh, nnzLocal: ", &
            no_l, numColLocal, maxnh, nnzLocal
  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo

else if (mpirank == mpisize-1) then

  numColLocal = bs + norbs_slack
  allocate(tmpi(norbs_slack))
  call MPI_Recv(tmpi,norbs_slack,MPI_Integer,0,0,MPI_comm,stat,ierr)
  call MPI_GET_COUNT(stat, MPI_Integer, count, mpierror)
  print *, "Last NODE received ", count, &
           " records in tmpi. Expected: ", norbs_slack

   nnzold = sum(numh(1:bs))
  if (nnzold /= maxnh) call die("nnz_old /= maxnh in last node")
  nnz_slack = sum(tmpi(1:norbs_slack))  
  nnzLocal = nnzold + nnz_slack
  print *, "Last node: no_l, numColLocal, maxnh, nnzLocal: ", no_l, numColLocal, maxnh, nnzLocal

  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,bs
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo
  do ih = bs+1,numColLocal
     colptrLocal(ih+1) = colptrLocal(ih) + tmpi(ih)
  enddo
  deallocate(tmpi)

else

  numColLocal = bs 
  nnzLocal = sum(numh(1:bs))
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
   allocate(FDMnzvalLocal(1:nnzLocal))

else if (mpirank == mpisize-1) then
   nnzold = sum(numh(1:bs))
   allocate(rowindLocal(1:nnzLocal))
   allocate(HnzvalLocal(1:nnzLocal))
   allocate(SnzvalLocal(1:nnzLocal))
   allocate(DMnzvalLocal(1:nnzLocal))
   allocate(EDMnzvalLocal(1:nnzLocal))
   allocate(FDMnzvalLocal(1:nnzLocal))

   rowindLocal(1:nnzold) =  listh(1:nnzold)
   HnzvalLocal(1:nnzold) =  H(1:nnzold,ispin)
   SnzvalLocal(1:nnzold) =  S(1:nnzold)
   call MPI_Recv(rowindLocal(nnzold+1),nnz_slack,MPI_Integer, &
                 0,1,MPI_comm,stat,ierr)
   call MPI_GET_COUNT(stat, MPI_Integer, count, mpierror)
   print *, "Last NODE received ", count, &
           " records from listh. Expected: ", nnz_slack
   call MPI_Recv(SnzvalLocal(nnzold+1),nnz_slack,MPI_Double_Precision, &
                 0,2,MPI_comm,stat,ierr)
   call MPI_Recv(HnzvalLocal(nnzold+1),nnz_slack,MPI_Double_Precision, &
                 0,3,MPI_comm,stat,ierr)

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


!  DMnzvalLocal(:) = 2*HnzvalLocal(:)
!  EDMnzvalLocal(:) = -1.0_dp

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
 call MPI_Recv(DM(nnzLocal+1,ispin),nnz_slack,MPI_Double_Precision, &
               mpisize-1,1,MPI_comm,stat,ierr)
 call MPI_Recv(EDM(nnzLocal+1,ispin),nnz_slack,MPI_Double_Precision, &
               mpisize-1,2,MPI_comm,stat,ierr)

 print *, "Node: ", mpirank, ": H: ",  H(nnzLocal+1:nnzLocal+2,ispin)
 print *, "Node: ", mpirank, ": DM: ", DM(nnzLocal+1:nnzLocal+2,ispin)
 print *, "Node: ", mpirank, ": EDM: ", EDM(nnzLocal+1:nnzLocal+2,ispin)

else if (mpirank == mpisize-1) then

   call MPI_Send(DMnzvalLocal(nnzold+1),nnz_slack,MPI_Double_Precision, &
                 0,1,MPI_comm,ierr)
   call MPI_Send(EDMnzvalLocal(nnzold+1),nnz_slack,MPI_Double_Precision, &
                 0,2,MPI_comm,ierr)
   DM(1:nnzold,ispin) = DMnzvalLocal(1:nnzold)
   EDM(1:nnzold,ispin) = EDMnzvalLocal(1:nnzold)

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

  end subroutine pexsi_EXCHANGE
end module m_pexsi_exchange
