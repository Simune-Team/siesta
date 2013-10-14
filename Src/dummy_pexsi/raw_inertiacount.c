#define FORTRAN(name) name##_
void FORTRAN(f_ppexsi_raw_inertiacount_interface)(
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
double*       eMin,                         // Lower bound for energy interval
double*       eMax,                         // Upper bound for energy interval
int*          nShifts,                      // Number of shifts in computing the inertia
int*          ordering,                     // SuperLUDIST ordering
int*          npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
int*    	    Fcomm,                        // Overall MPI communicator
// Output parameters
double*       muList,                       // The list of shifts
int*       numElectronList,              // The number of electrons (raw integer counts) corresponding to shifts
int*          info                          // 0: successful exit.  1: unsuccessful
		)
{

	return;
} // -----  end of function f_ppexsi_raw_inertiacount_interface  ----- 

