module m_pexsi
  public :: pexsi
CONTAINS
  subroutine pexsi(no_u, nspin,  &
       maxnh, numh, listhptr, listh, H, S, qtot)

    use precision, only  : dp
    use m_hsx, only      : write_hs_formatted
    use m_walltime, only : wall_time
#ifdef MPI
    use mpi_siesta, only: mpi_comm_world
#endif
    implicit          none

    integer, intent(in)  :: maxnh, no_u, nspin
    integer, intent(in)  :: listh(maxnh), numh(*), listhptr(*)
    real(dp), intent(in) :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: qtot

!Lin variables
integer :: nrows, nnz, nnzLocal, numColLocal
integer, allocatable, dimension(:) ::  colptrLocal, rowindLocal
real(dp), allocatable, dimension(:) :: &
	HnzvalLocal, SnzvalLocal, DMnzvalLocal, EDMnzvalLocal, &
	FDMnzvalLocal
integer :: numPole
real(dp) :: temperature, numElectronExact, numElectron,&
	gap, deltaE
real(dp) :: mu, muMin, muMax
integer:: muMaxIter
real(dp) :: poleTolerance, numElectronTolerance
integer:: npPerPole, nprow, npcol
integer :: mpirank, mpisize, ierr
real(dp):: timeSta, timeEnd
character*32 :: Hfile, Sfile
!------------


    external      :: timer

    call timer("pexsi", 1)

    call write_hs_formatted(no_u, nspin,  &
         maxnh, numh, listhptr, listh, H, S)


!!!!call mpi_init( ierr )
call mpi_comm_rank( MPI_COMM_WORLD, mpirank, ierr )
call mpi_comm_size( MPI_COMM_WORLD, mpisize, ierr )


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
Hfile            = "H.matrix"
Sfile            = "S.matrix"

! Read and compute the size/local size of the arrays 
! The conversion of the string to fit the C format is important.
call f_read_distsparsematrix_formatted_head( &
	trim(Hfile)//char(0),&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	MPI_COMM_WORLD )

if( mpirank .eq. 0 ) then
	write(*,*) "Matrix size (local data on proc 0):" 
	write(*,*) "size = ", nrows
	write(*,*) "nnz  = ", nnz
	write(*,*) "nnzLocal = ", nnzLocal
	write(*,*) "numColLocal = ", numColLocal
endif

! Allocate memory
allocate( colptrLocal( numColLocal + 1 ) )
allocate( rowindLocal( nnzLocal ) )
allocate( HnzvalLocal( nnzLocal ) )
allocate( SnzvalLocal( nnzLocal ) ) 
allocate( DMnzvalLocal( nnzLocal ) ) 
allocate( EDMnzvalLocal( nnzLocal ) ) 
allocate( FDMnzvalLocal( nnzLocal ) ) 

!timeSta = mpi_wtime()
call wall_time(timeSta)

call f_read_distsparsematrix_formatted (&
	trim(Hfile)//char(0),&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
	rowindLocal,&
	HnzvalLocal,&
	MPI_COMM_WORLD )

call f_read_distsparsematrix_formatted (&
	trim(Sfile)//char(0),&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
	rowindLocal,&
	SnzvalLocal,&
	MPI_COMM_WORLD )

!timeEnd = mpi_wtime()
call wall_time(timeEnd)

if( mpirank == 0 ) then
  write(*,*) "Time for reading H/S matrices is ", &
		timeEnd - timeSta, " [s]"
endif

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
	poleTolerance,&
	numElectronTolerance,&
	MPI_COMM_WORLD,&
	npPerPole )

if( mpirank == 0 ) then
	write(*, *) "mu          = ", mu
	write(*, *) "numElectron = ", numElectron
endif


deallocate( colptrLocal )
deallocate( rowindLocal )
deallocate( HnzvalLocal )
deallocate( SnzvalLocal )
deallocate( DMnzvalLocal )
deallocate( EDMnzvalLocal )
deallocate( FDMnzvalLocal )

    call timer("pexsi", 2)

  end subroutine pexsi
end module m_pexsi
