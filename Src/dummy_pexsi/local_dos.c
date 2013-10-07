#define FORTRAN(name) name##_

/// @brief Interface between PPEXSI and C for computing the local
/// density of states, i.e.
///   n(r;E) = lim_{eta->0+} Im 1/pi <r| (H-(E+i eta)I)^{-1} |r>
/// This code returns the selected elements of the matrix
///   
///   1/pi Im ( H - (E+i eta) S )^{-1}
///   
void FORTRAN(f_ppexsi_localdos_interface)(
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
		double*       Energy,                       // Real part of the shift
		double*       eta,                          // Broadening parameter
		int*          ordering,                     // SuperLUDIST ordering
		int*          npSymbFact,                   // Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
		int*    	    Fcomm,                        // Overall MPI communicator
		// Output parameters
		double*       localDOSnzvalLocal,           // Nonzero value of Im 1/pi (H - (E+ieta) S)^{-1}
		int*          info                          // 0: successful exit.  1: unsuccessful
		)
{
  return;
} // -----  end of function f_ppexsi_localdos_interface  ----- 
