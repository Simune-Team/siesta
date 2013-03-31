!
! Header for the interface block for PEXSI
!
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

