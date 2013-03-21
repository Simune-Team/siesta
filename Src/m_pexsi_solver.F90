module m_pexsi_solver

  public :: pexsi_solver

CONTAINS

! This version uses the standard PEXSI distribution (by tricking Siesta into
! using it)
!
  subroutine pexsi_solver(no_u, no_l, nspin,  &
       maxnh, numh, listhptr, listh, H, S, qtot, DM, EDM, &
       ef, freeEnergyCorrection, temp)

    use precision, only  : dp
    use fdf
    use parallel, only   : worker, ionode
    use m_mpi_utils, only: globalize_sum
    use units,       only: Kelvin, eV
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

    integer  :: ordering
    integer  :: muIter
    real(dp) :: muZeroT

    real(dp), save :: mu
    logical, save  :: first_call = .true.
    real(dp)       :: eBandStructure, eBandH

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
real(dp), allocatable, dimension(:) :: muList, &
                                       numElectronList, &
                                       numElectronDrvList
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

interface
 subroutine f_ppexsi_interface( &
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
        rowindLocal,&
	HnzvalLocal,&
	SnzvalLocal,&
	DMnzvalLocal,&
	EDMnzvalLocal,&
	FDMnzvalLocal,&
	numPole,&
	temperature,&
	numElectronExact,&
	numElectron,&
	gap,&
	deltaE,&
	mu,&
	muMin,&
        muMax,&
	muMaxIter,&
        ordering, &
        muIter, &
        muList, &
        numElectronList,&
        numElectronDrvList,&
        muZeroT,&
	poleTolerance,&
	numElectronTolerance,&
	comm_global,&
	npPerPole )

   integer, intent(in) :: nrows, nnz, nnzLocal, numColLocal
   integer, intent(in) :: colptrLocal(:), rowindLocal(:)
   real(SELECTED_REAL_KIND(10,100)), intent(in) :: HnzvalLocal(:),&
                                                   SnzvalLocal(:)
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: DMnzvalLocal(:),&
                                                    EDMnzvalLocal(:), &
                                                    FDMnzvalLocal(:)
   integer, intent(in)                            :: numPole
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: temperature
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: numElectronExact
   real(SELECTED_REAL_KIND(10,100)), intent(out)  :: numElectron
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: gap, deltaE
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: muMin, muMax
   real(SELECTED_REAL_KIND(10,100)), intent(inout):: mu

   ! Variables related to mu history
   ! Maximum number of allowed iterations
   integer, intent(in)                           :: muMaxIter

   ! Ordering 
   !   0   : PARMETIS
   !   1   : METIS_AT_PLUS_A
   !   2   : MMD_AT_PLUS_A
   integer, intent(in)                           :: ordering

   ! Actual number of iterations performed
   integer, intent(out)                          :: muIter

   ! List of values of mu, N_e, d(N_e)/d_mu
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: muList(muMaxIter)
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: numElectronList(muMaxIter)
   real(SELECTED_REAL_KIND(10,100)),intent(out):: numElectronDrvList(muMaxIter)

   ! mu extrapolated to T=0K
   real(SELECTED_REAL_KIND(10,100)), intent(out) :: muZeroT


   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: poleTolerance, &
                                                     numElectronTolerance
   integer, intent(in)                 :: comm_global, npPerPole

 end subroutine f_ppexsi_interface
end interface

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
   temperature      = temp/Kelvin
   if (IOnode) write(6,*) "Electronic temperature: ", temperature

   call MPI_Barrier(Siesta_comm,ierr)

   nrows = no_u
!  numColLocal = num_local_elements(dist2,no_u,myrank2)
   numColLocal = no_l

   ! Figure out the communication needs
   ! call analyze_comms(comms)
   ! call do_transfers(comms,numh,numh_pexsi, &
   !                   g1,g2,mpi_comm)

 !   nnzLocal = sum(numh_pexsi(1:numColLocal))
 ! call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,PEXSI_comm,MPIerror)
   nnzLocal = sum(numh(1:numColLocal))
   call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,Siesta_comm,MPIerror)


!  print *, "Node ", mpirank, ": no_l, numColLocal, maxnh, nnzLocal: ", &
!           no_l, numColLocal, maxnh, nnzLocal

  allocate(colptrLocal(1:numColLocal+1))
  colptrLocal(1) = 1
  do ih = 1,numColLocal
!    colptrLocal(ih+1) = colptrLocal(ih) + numh_pexsi(ih)
     colptrLocal(ih+1) = colptrLocal(ih) + numh(ih)
  enddo

   allocate(rowindLocal(1:nnzLocal))
   allocate(HnzvalLocal(1:nnzLocal))
   allocate(SnzvalLocal(1:nnzLocal))
   allocate(DMnzvalLocal(1:nnzLocal))
   allocate(EDMnzvalLocal(1:nnzLocal))
   allocate(FDMnzvalLocal(1:nnzLocal))

!  call generate_commsnnz(comms,commsnnz)
!  call do_transfers(commsnnz,listh,rowindLocal, &
!                        g1, g2, mpi_comm)
!  call do_transfers(commsnnz,H,HnzvalLocal, &
!                        g1, g2, mpi_comm)
!  call do_transfers(commsnnz,S,SnzvalLocal, &
!                        g1, g2, mpi_comm)

   rowindLocal(1:nnzLocal) = listh(1:nnzLocal)
   HnzvalLocal(1:nnzLocal) = H(1:nnzLocal,ispin)
   SnzvalLocal(1:nnzLocal) = S(1:nnzLocal)

endif ! worker

!temperature      = fdf_get("PEXSI.temperature",3000.0d0)    ! Units??
! Now passed directly by Siesta  (Use ElectronicTemperature (with units))

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

allocate( muList( muMaxIter ) )
allocate( numElectronList( muMaxIter ) )
allocate( numElectronDrvList( muMaxIter ) )

! Do not compute a pole if the corresponding weight is < poleTolerance.
poleTolerance    = fdf_get("PEXSI.pole-tolerance",1d-8)

! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
numElectronTolerance = fdf_get("PEXSI.num-electron-tolerance",1d-1)

! Number of processors used for each pole. At the moment use mpisize.
! Later can be changed to 
!npPerPole        = fdf_get("PEXSI.np-per-pole",mpisize)

! Ordering flag
ordering = fdf_get("PEXSI.ordering",1)

call MPI_Bcast(npPerPole,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(nrows,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(numElectronExact,1,MPI_double_precision,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(temperature,1,MPI_double_precision,0,true_MPI_COMM_world,ierr)
call MPI_Bcast(ordering,1,MPI_integer,0,true_MPI_COMM_world,ierr)

call f_ppexsi_interface( &
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
        rowindLocal,&
	HnzvalLocal,&
	SnzvalLocal,&
	DMnzvalLocal,&
	EDMnzvalLocal,&
	FDMnzvalLocal,&
	numPole,&
	temperature,&
	numElectronExact,&
	numElectron,&
	gap,&
	deltaE,&
	mu,&
	muMin,&
        muMax,&
	muMaxIter,&
        ordering, &
        muIter, &
        muList, &
        numElectronList,&
        numElectronDrvList,&
        muZeroT,&
	poleTolerance,&
	numElectronTolerance,&
	true_MPI_COMM_WORLD,&
	npPerPole )

!  call reverse_comms(commsnnz,commsnnz_reverse)
!  call do_transfers(commsnnz_reverse,DMnzvalLocal,DM, &
!                        g2, g1, mpi_comm)
!  call do_transfers(commsnnz_reverse,EDMnzvalLocal,EDM, &
!                        g2, g1, mpi_comm)
!  allocate(FDM(1:maxnh))
!  call do_transfers(commsnnz_reverse,FDMnzvalLocal,FDM, &
!                        g2, g1, mpi_comm)

if (worker) then

   ! Watch out here. We have to assume that rank 0 is common to both
   ! subgroups
   ! if (mpirank == 0) then
   !  ef = mu
   !  call MPI_Bcast(ef, .... siesta subgroup)
   ! endif

   ef = mu

   ! this will disappear
   DM(1:nnzLocal,ispin) = DMnzvalLocal(1:nnzLocal)
   EDM(1:nnzLocal,ispin) = EDMnzvalLocal(1:nnzLocal)

   freeEnergyCorrection = 0.0_dp
   eBandStructure = 0.0_dp
   eBandH = 0.0_dp
   do i = 1,nnzLocal
      ! Watch out for FDM nzvalLocal (see above)
      freeEnergyCorrection = freeEnergyCorrection + SnzvalLocal(i) * &
           ( FDMnzvalLocal(i) - EDMnzvalLocal(i) )
      eBandStructure = eBandStructure + SnzvalLocal(i) * &
           ( EDMnzvalLocal(i) )
      eBandH = eBandH + HnzvalLocal(i) * &
           ( DMnzvalLocal(i) )
   enddo

   ! This operations in siesta group (implicit communicator)
   call globalize_sum( freeEnergyCorrection, buffer1 )
   freeEnergyCorrection = buffer1 + mu*numElectron
   call globalize_sum( eBandStructure, buffer1 )
   eBandStructure = buffer1
   call globalize_sum( eBandH, buffer1 )
   eBandH = buffer1

   if( mpirank == 0 ) then
      write(*, *) "mu          = ", mu
      write(*, *) "mu (eV)     = ", mu/eV
      write(*, *) "muZeroT (eV)     = ", muZeroT/eV
      write(*, *) "numElectron = ", numElectron
      write(*, *) "eBandStructure (Ry) = ", eBandStructure
      write(*, *) "eBandStructure (eV) = ", eBandStructure/eV
      write(*, *) "eBandH (eV) = ", eBandH/eV
      write(*, *) "freeEnergy (eV) = ", (eBandStructure + freeEnergyCorrection)/eV
      write(*,*) "Number of mu iterations: ", muIter
      write(*,"(a3,2a12,a20)") "it", "mu", "N_e", "dN_e/dmu"
      do i = 1, muIter
         write(*,"(i3,2f12.4,g20.5)") i, muList(i)/eV, &
                 numElectronList(i), numElectronDrvList(i)*eV
      end do
   endif

  deallocate(rowindLocal)
  deallocate(HnzvalLocal)
  deallocate(SnzvalLocal)
  deallocate(DMnzvalLocal)
  deallocate(EDMnzvalLocal)
  deallocate(FDMnzvalLocal)
  deallocate(colPtrLocal)
! deallocate(FDM)

    call timer("pexsi", 2)
 endif ! worker
#endif 

end subroutine pexsi_solver
end module m_pexsi_solver
