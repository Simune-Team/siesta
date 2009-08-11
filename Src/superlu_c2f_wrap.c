#ifdef _LANCZOS_
#include "superlu_ddefs.h"
#include "Cnames.h"
#endif
/* functions that create memory for a struct and return a handle */

#if 0
typedef int fptr;  /* 32-bit */
#else
typedef long int fptr;  /* 64-bit */
#endif


#ifdef _LANCZOS_
/* functions that create memory for a struct and return a handle */
void f_create_gridinfo_handle(fptr *handle)
{
  *handle = (fptr) SUPERLU_MALLOC(sizeof(gridinfo_t));
}

void f_create_options_handle(fptr *handle)
{
  *handle = (fptr) SUPERLU_MALLOC(sizeof(superlu_options_t));
}

void f_create_ScalePerm_handle(fptr *handle)
{
  *handle = (fptr) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
}

void f_create_LUstruct_handle(fptr *handle)
{
  *handle = (fptr) SUPERLU_MALLOC(sizeof(LUstruct_t));
}

void f_create_SOLVEstruct_handle(fptr *handle)
{
  *handle = (fptr) SUPERLU_MALLOC(sizeof(SOLVEstruct_t));
}

void f_create_SuperMatrix_handle(fptr *handle)
{
  *handle = (fptr) SUPERLU_MALLOC(sizeof(SuperMatrix));
}

void f_create_SuperLUStat_handle(fptr *handle)
{
  *handle = (fptr) SUPERLU_MALLOC(sizeof(SuperLUStat_t));
}

/* functions that free the memory allocated by the above functions */
void f_destroy_gridinfo_handle(fptr *handle)
{
  SUPERLU_FREE((void *)*handle);
}

void f_destroy_options_handle(fptr *handle)
{
  SUPERLU_FREE((void *)*handle);
}

void f_destroy_ScalePerm_handle(fptr *handle)
{
  SUPERLU_FREE((void *)*handle);
}

void f_destroy_LUstruct_handle(fptr *handle)
{
  SUPERLU_FREE((void *)*handle);
}

void f_destroy_SOLVEstruct_handle(fptr *handle)
{
  SUPERLU_FREE((void *)*handle);
}

void f_destroy_SuperMatrix_handle(fptr *handle)
{
  SUPERLU_FREE((void *)*handle);
}

void f_destroy_SuperLUStat_handle(fptr *handle)
{
  SUPERLU_FREE((void *)*handle);
}

/* functions that get or set values in a C struct.
   This is not the complete set of structs for which a user might want
   to get/set a component, and there may be missing components. */

void f_get_gridinfo(fptr *grid, int *iam, int *nprow, int *npcol)
{
  *iam=((gridinfo_t *) *grid)->iam;
  *npcol=((gridinfo_t *) *grid)->npcol;
  *nprow=((gridinfo_t *) *grid)->nprow;
}

void f_get_SuperMatrix(fptr *A, int *nrow, int *ncol)
{
  *nrow = ((SuperMatrix *) *A)->nrow;
  *ncol = ((SuperMatrix *) *A)->ncol;
}

void f_set_SuperMatrix(fptr *A, int *nrow, int *ncol)
{
  ((SuperMatrix *) *A)->nrow = *nrow;
  ((SuperMatrix *) *A)->ncol = *ncol;
}

void f_get_CompRowLoc_Matrix(fptr *A, int *m, int *n, int *nnz_loc,
			     int *m_loc, int *fst_row)
{
  *m=((SuperMatrix *) *A)->nrow;
  *n=((SuperMatrix *) *A)->ncol;
  *m_loc=((NRformat_loc *) ((SuperMatrix *) *A)->Store)->m_loc;
  *nnz_loc=((NRformat_loc *) ((SuperMatrix *) *A)->Store)->nnz_loc;
  *fst_row=((NRformat_loc *) ((SuperMatrix *) *A)->Store)->fst_row;
}

void f_set_CompRowLoc_Matrix(fptr *A, int *m, int *n, int *nnz_loc,
			     int *m_loc, int *fst_row)
{
  ((SuperMatrix *) *A)->nrow = *m;
  ((SuperMatrix *) *A)->ncol = *n;
  ((NRformat_loc *) ((SuperMatrix *) *A)->Store)->m_loc = *m_loc;
  ((NRformat_loc *) ((SuperMatrix *) *A)->Store)->nnz_loc = *nnz_loc;
  ((NRformat_loc *) ((SuperMatrix *) *A)->Store)->fst_row = *fst_row;
}

void f_get_superlu_options(fptr *opt, int *Fact, int *Equil, int *ParSymbFact,
                           int *ColPerm, int *RowPerm, int *IterRefine,
			   int *Trans, int *ReplaceTinyPivot,
			   int *SolveInitialized, int *RefineInitialized,
			   int *PrintStat)
{
   *Fact = (int) ((superlu_options_t *) *opt)->Fact;
   *Equil = (int) ((superlu_options_t *) *opt)->Equil;
   *ParSymbFact = (int) ((superlu_options_t *) *opt)->ParSymbFact;
   *ColPerm = (int) ((superlu_options_t *) *opt)->ColPerm;
   *RowPerm = (int) ((superlu_options_t *) *opt)->RowPerm;
   *IterRefine = (int) ((superlu_options_t *) *opt)->IterRefine;
   *Trans = (int) ((superlu_options_t *) *opt)->Trans;
   *ReplaceTinyPivot = (int) ((superlu_options_t *) *opt)->ReplaceTinyPivot;
   *SolveInitialized = (int) ((superlu_options_t *) *opt)->SolveInitialized;
   *RefineInitialized = (int) ((superlu_options_t *) *opt)->RefineInitialized;
   *PrintStat = (int) ((superlu_options_t *) *opt)->PrintStat;
}

void f_set_superlu_options(fptr *opt, int *Fact, int *Equil, int *ParSymbFact,
                           int *ColPerm, int *RowPerm, int *IterRefine,
			   int *Trans, int *ReplaceTinyPivot,
			   int *SolveInitialized, int *RefineInitialized,
			   int *PrintStat)
{
    superlu_options_t *l_options = (superlu_options_t*) *opt;
    l_options->Fact = (fact_t) *Fact;
   ((superlu_options_t *) *opt)->Equil = (yes_no_t) *Equil;
   ((superlu_options_t *) *opt)->ParSymbFact = (yes_no_t) *ParSymbFact;
   ((superlu_options_t *) *opt)->ColPerm = (colperm_t) *ColPerm;
   ((superlu_options_t *) *opt)->RowPerm = (rowperm_t) *RowPerm;
   ((superlu_options_t *) *opt)->IterRefine = (IterRefine_t) *IterRefine;
   ((superlu_options_t *) *opt)->Trans = (trans_t) *Trans;
   ((superlu_options_t *) *opt)->ReplaceTinyPivot = (yes_no_t) *ReplaceTinyPivot;
   ((superlu_options_t *) *opt)->SolveInitialized = (yes_no_t) *SolveInitialized;
   ((superlu_options_t *) *opt)->RefineInitialized = (yes_no_t) *RefineInitialized;
   ((superlu_options_t *) *opt)->PrintStat = (yes_no_t) *PrintStat;
}

/* wrappers for SuperLU functions */

void f_set_default_options(fptr *options)
{
   set_default_options_dist((superlu_options_t *) *options);
}

void f_superlu_gridinit(int *Bcomm, int *nprow, int *npcol, fptr *grid)
{
   superlu_gridinit( (int_t) *Bcomm, (int_t) *nprow, (int_t) *npcol,
                    (gridinfo_t *) *grid);
}

void f_superlu_gridexit(fptr *grid)
{
   superlu_gridexit((gridinfo_t *) *grid);
}

void f_ScalePermstructInit(int *m, int *n, fptr *ScalePermstruct)
{
   ScalePermstructInit((int_t) *m, (int_t) *n,
                       (ScalePermstruct_t *) *ScalePermstruct);
}

void f_ScalePermstructFree(fptr *ScalePermstruct)
{
   ScalePermstructFree((ScalePermstruct_t *) *ScalePermstruct);
}

void f_PStatInit(fptr *stat)
{
   PStatInit((SuperLUStat_t *) *stat);
}

void f_PStatFree(fptr *stat)
{
   PStatFree((SuperLUStat_t *) *stat);
}

void f_LUstructInit(int *m, int *n, fptr *LUstruct)
{
   LUstructInit((int_t) *m, (int_t) *n, (LUstruct_t *) *LUstruct);
}

void f_LUstructFree(fptr *LUstruct)
{
   LUstructFree((LUstruct_t *) *LUstruct);
}

void f_Destroy_LU(int *n, fptr *grid, fptr *LUstruct)
{
   Destroy_LU((int_t) *n, (gridinfo_t *) *grid, (LUstruct_t *) *LUstruct);
}

void f_dCreate_CompRowLoc_Mat_dist(fptr *A, int *m, int *n, int *nnz_loc,
				   int *m_loc, int *fst_row, double *nzval,
				   int *colind, int *rowptr, int *stype,
				   int *dtype, int *mtype)
{
   dCreate_CompRowLoc_Matrix_dist((SuperMatrix *) *A, (int_t) *m, (int_t) *n,
                                  (int_t) *nnz_loc, (int_t) *m_loc,
                                  (int_t) *fst_row, (double *) nzval,
                                  (int_t *) colind, (int_t *) rowptr,
                                  (Stype_t) *stype, (Dtype_t) *dtype,
                                  (Mtype_t) *mtype);
}

void f_Destroy_CompRowLoc_Mat_dist(fptr *A)
{
   Destroy_CompRowLoc_Matrix_dist((SuperMatrix *) *A);
}

void f_Destroy_SuperMat_Store_dist(fptr *A)
{
   Destroy_SuperMatrix_Store_dist((SuperMatrix *) *A);
}

void f_dSolveFinalize(fptr *options, fptr *SOLVEstruct)
{
   dSolveFinalize((superlu_options_t *) *options,
                  (SOLVEstruct_t *) *SOLVEstruct);
}

void f_pdgssvx(fptr *options, fptr *A, fptr *ScalePermstruct, double *B,
               int *ldb, int *nrhs, fptr *grid, fptr *LUstruct,
               fptr *SOLVEstruct, double *berr, fptr *stat, int *info)
{
    pdgssvx((superlu_options_t *) *options, (SuperMatrix *) *A,
	    (ScalePermstruct_t *) *ScalePermstruct, B, *ldb, *nrhs,
	    (gridinfo_t *) *grid, (LUstruct_t *) *LUstruct,
	    (SOLVEstruct_t *) *SOLVEstruct, berr,
	    (SuperLUStat_t *) *stat, info);

    PStatPrint((superlu_options_t *) *options, (SuperLUStat_t *) *stat,
	       (gridinfo_t *) *grid);
}
#endif
