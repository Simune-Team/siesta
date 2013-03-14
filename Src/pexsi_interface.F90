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
   integer, intent(in)                 :: numPole
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: temperature
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: numElectronExact
   real(SELECTED_REAL_KIND(10,100)), intent(out)  :: numElectron
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: gap, deltaE, muMin, muMax
   real(SELECTED_REAL_KIND(10,100)), intent(inout):: mu
   integer, intent(in)                            :: muMaxIter
   real(SELECTED_REAL_KIND(10,100)), intent(in)   :: poleTolerance, &
                                                     numElectronTolerance
   integer, intent(in)                 :: comm_global, npPerPole

 end subroutine f_ppexsi_interface
