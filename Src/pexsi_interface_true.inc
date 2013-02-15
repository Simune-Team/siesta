
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
	true_MPI_COMM_WORLD,&
	npPerPole )
