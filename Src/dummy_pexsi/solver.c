#define FORTRAN(name) name##_

/// @brief Interface between PPEXSI and FORTRAN for the solve procedure.
void FORTRAN(f_ppexsi_solve_interface)(
		// Input parameters
		int*          nrows,                        // Size of the matrix
	  int*          nnz,                          // Total number of nonzeros in H
		int*          nnzLocal,                     // Number of nonzeros in H on this proc
		int*          numColLocal,                  // Number of local columns for H
		int*          colptrLocal,                  // Colomn pointer in CSC format
		int*          rowindLocal,                  // Row index pointer in CSC format
		double*       HnzvalLocal,                  // Nonzero value of H in CSC format
		int*          isSIdentity,                  // Whether S is an identity matrix. If so, the variable SnzvalLocal is omitted.
		double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
		double*       temperature,                  // Temperature, in the same unit as H
		double*       numElectronExact,             // Exact number of electrons
		double*       mu0,                          // Initial guess for mu
		double*       muMin0,                       // Initial guess for lower bound of mu
		double*       muMax0,                       // Initial guess for upper bound of mu
		double*       gap,                          // Energy gap (lower bound)
		double*       deltaE,                       // Spectral radius of S^{-1}H
		int*          numPole,                      // Number of poles
		int*          maxIter,                      // Maximum number of iterations for mu-iteration in PEXSI
		double*       numElectronTolerance,         // Stopping criterion of PEXSI mu iteration.
		int*          ordering,                     // SuperLUDIST ordering
	  int*          npPerPole,                    // Number of processors for each pole
		int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
		int*    	    Fcomm,                        // Overall MPI communicator
		// Output parameters
		double*      DMnzvalLocal,                  // Nonzero value of density matrix in CSC format
		double*     EDMnzvalLocal,                  // Nonzero value of energy density matrix in CSC format
		double*     FDMnzvalLocal,                  // Nonzero value of free energy density matrix in CSC format
		double*       muPEXSI,                      // Final chemical potential
		double*       numElectronPEXSI,             // Computed number of electron at the final chemical potential
		double*       muMinPEXSI,                   // Final lower bound for mu.
		double*       muMaxPEXSI,                   // Final upper bound for mu
		int*          numIter,                      // Number of actual iterations for PEXSI
		double*       muList,                       // The history of mu
		double*       numElectronList,              // The history of number of electrons correspondig to mu
		double*       numElectronDrvList,           // The history of dN/dMu
		int*          info                          // 0: successful exit.  1: unsuccessful
		)
{
	return;
} // -----  end of function f_ppexsi_solve_interface  ----- 

