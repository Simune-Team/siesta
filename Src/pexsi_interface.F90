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


   integer, parameter :: dp = SELECTED_REAL_KIND(10,100)

   integer :: ierr, mpirank

   call MPI_Comm_Rank(comm_global, mpirank, ierr)

   if (mpirank == 0) then
      write(*,*) "PEXSI:"
      write(*,*) "nrows: ", nrows
      write(*,*) "nnz: ", nnz
      write(*,*) "nnzLocal (0): ", nnzLocal
      write(*,*) "numcolLocal (0): ", numColLocal
      write(*,*) "numPole: ", numPole
      write(*,*) "temperature: ", temperature
      write(*,*) "numElectronExact: ", numElectronExact
      write(*,*) "numElectronTolerance: ", numElectronTolerance
      write(*,*) "gap: ", gap
      write(*,*) "DeltaE: ", deltaE
      write(*,*) "muMin: ", muMin
      write(*,*) "muMax: ", muMax
      write(*,*) "muMaxIter: ", muMaxIter
      write(*,*) "npPerPole: ", npPerPole
      write(*,*) "PoleTolerance: ", poleTolerance
      write(*,*) "ordering: ", ordering
   endif

!  Sample dummy outputs

   muIter = min(muMaxIter,2)

   numElectron = 24.0_dp
   mu = -0.3_dp

   numElectronList(1:muIter-1) = 23.9_dp
   numElectronList(muIter) = numElectron

   muList(1:muIter-1) = -0.29_dp
   muList(muIter) = mu

   numElectronDrvList(1:muIter) = 0.04_dp

   muZeroT = -0.28_dp
   

 end subroutine f_ppexsi_interface
