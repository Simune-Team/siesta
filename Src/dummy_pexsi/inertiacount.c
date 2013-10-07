#define FORTRAN(name) name##_
void FORTRAN(f_ppexsi_inertiacount_interface)(
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
double*       muMin0,                       // Initial guess of lower bound for mu
double*       muMax0,                       // Initial guess of upper bound for mu
int*          numPole,                      // Number of shifts in computing the inertia, still called "Pole" for legacy reason
int*          maxIter,                      // Maximum number of iterations for computing the inertia
double*       numElectronTolerance,         // Stopping criterion of inertia count
int*          ordering,                     // SuperLUDIST ordering
int*          npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
int*    	    Fcomm,                        // Overall MPI communicator
// Output parameters
double*       muMinInertia,                 // Lower bound for mu after inertia count
double*       muMaxInertia,                 // Upper bound for mu after inertia count
double*       muLowerEdge,                  // Ne(muLowerEdge) = Ne - eps. For band gapped system
double*       muUpperEdge,                  // Ne(muUpperEdge) = Ne + eps. For band gapped system
int*          numIter,                      // Number of actual iterations for inertia count
double*       muList,                       // The list of shifts
double*       numElectronList,              // The number of electrons (finite temperature) corresponding to shifts
int*          info                          // 0: successful exit.  1: unsuccessful
		)
{

	return;
} // -----  end of function f_ppexsi_inertiacount_interface  ----- 

