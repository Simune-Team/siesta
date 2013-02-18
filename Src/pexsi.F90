module m_pexsi

  public :: pexsi

CONTAINS

! SIESTA interface to PEXSI
!
! This version uses the standard PEXSI distribution (by tricking Siesta into
! using it)
!
  subroutine pexsi(no_u, no_l, nspin,  &
       maxnh, numh, listhptr, listh, H, S, qtot, DM, EDM, &
       ef, freeEnergyCorrection, temp)

    use precision, only  : dp
    use fdf
    use parallel, only   : worker
    use m_mpi_utils, only: globalize_sum
    use units,       only: Kelvin
#ifdef MPI
    use mpi_siesta
#endif
    implicit          none

    integer, intent(in)  :: maxnh, no_u, no_l, nspin
    integer, intent(in), target  :: listh(maxnh), numh(*), listhptr(*)
    real(dp), intent(in), target :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: qtot
    real(dp), intent(out), target:: DM(maxnh,nspin), EDM(maxnh,nspin)
    real(dp), intent(out)        :: ef  ! Fermi energy
    real(dp), intent(out)        :: freeEnergyCorrection
    real(dp), intent(in)         :: temp   ! Electronic temperature

#ifndef MPI
    call die("PEXSI needs MPI")
#else

    integer :: siesta_comm 
    integer :: MPIerror, stat(MPI_STATUS_SIZE), count
    integer :: bs, norbs_slack, nnz_slack
    integer :: ispin, maxnhtot, ih, nnzold, i

    real(dp), save :: mu
    logical, save  :: first_call = .false.

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
real(dp) :: muMin, muMax
integer:: muMaxIter
real(dp) :: poleTolerance, numElectronTolerance
integer:: npPerPole, nprow, npcol
integer :: mpirank, mpisize, ierr
!------------

real(dp) :: buffer1

external      :: timer

if (worker) then
   siesta_comm = mpi_comm_world
   call timer("pexsi", 1)

   ispin = 1
   if (nspin /=1) then
      call die("Spin polarization not yet supported in PEXSI")
   endif

   call mpi_comm_rank( SIESTA_COMM, mpirank, ierr )
   call mpi_comm_size( SIESTA_COMM, mpisize, ierr )

   npPerPole = mpisize
   numElectronExact = qtot   ! 2442.0d0 for DNA

   call MPI_Barrier(Siesta_comm,ierr)

   nrows = no_u
   numColLocal = no_l

   nnzLocal = sum(numh(1:numColLocal))
   call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,Siesta_comm,MPIerror)


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

endif ! worker

! Data is for the DNA matrix.

!temperature      = fdf_get("PEXSI.temperature",3000.0d0)    ! Units??
! Now passed directly by Siesta  (Use ElectronicTemperature (with units))
temperature      = temp/Kelvin
numPole          = fdf_get("PEXSI.num-poles",20)
gap              = fdf_get("PEXSI.gap",0.0d0)

! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
deltaE           = fdf_get("PEXSI.delta-E",3.0d0)

! Initial guess of chemical potential, also updated after pexsi.
if (first_call) then
   mu = fdf_get("PEXSI.mu",-0.60_dp)
   first_call = .false.
endif

! Lower/Upper bound for the chemical potential.
muMin            = fdf_get("PEXSI.mu-min",-1.0d0)
muMax            = fdf_get("PEXSI.mu-max", 0.0d0)

! muMaxIter should be 1 or 2 later when combined with SCF.
muMaxIter        = fdf_get("PEXSI.mu-max-iter",10)

! Do not compute a pole if the corresponding weight is < poleTolerance.
poleTolerance    = fdf_get("PEXSI.pole-tolerance",1d-8)

! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
numElectronTolerance = fdf_get("PEXSI.num-electron-tolerance",1d-1)

! Number of processors used for each pole. At the moment use mpisize.
! Later can be changed to 
!npPerPole        = fdf_get("PEXSI.np-per-pole",mpisize)

call MPI_Bcast(npPerPole,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(nrows,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(numElectronExact,1,MPI_double_precision,0,true_MPI_COMM_world,ierr)

include "pexsi_interface.inc"

if (worker) then
   if( mpirank == 0 ) then
      write(*, *) "mu          = ", mu
      write(*, *) "numElectron = ", numElectron
   endif

   ef = mu
   DM(1:nnzLocal,ispin) = DMnzvalLocal(1:nnzLocal)
   EDM(1:nnzLocal,ispin) = EDMnzvalLocal(1:nnzLocal)

   freeEnergyCorrection = 0.0_dp
   do i = 1,nnzLocal
      freeEnergyCorrection = freeEnergyCorrection + SnzvalLocal(i) * &
           ( FDMnzvalLocal(i) - EDMnzvalLocal(i) )
   enddo
   call globalize_sum( freeEnergyCorrection, buffer1 )
   freeEnergyCorrection = buffer1

  deallocate(rowindLocal)
  deallocate(HnzvalLocal)
  deallocate(SnzvalLocal)
  deallocate(DMnzvalLocal)
  deallocate(EDMnzvalLocal)
  deallocate(FDMnzvalLocal)
  deallocate(colPtrLocal)

    call timer("pexsi", 2)
 endif ! worker
#endif 

  end subroutine pexsi
end module m_pexsi
