      module SUPERLU
      use mpi
      implicit none

      integer, parameter :: cptr = 8 ! 64-bit
      integer, parameter :: NO                      = 0
      integer, parameter :: YES                     = 1
      integer, parameter :: DOFACT                  = 0 ! fact_t
      integer, parameter :: SamePattern             = 1
      integer, parameter :: SamePattern_SameRowPerm = 2
      integer, parameter :: FACTORED                = 3
      integer, parameter :: NOROWPERM               = 0 ! rowperm_t
      integer, parameter :: LargeDiag               = 1
      integer, parameter :: MY_PERMR                = 2
      integer, parameter :: NATURAL                 = 0 ! colperm_t
      integer, parameter :: MMD_AT_PLUS_A           = 1
      integer, parameter :: MMD_ATA                 = 2
      integer, parameter :: METIS_AT_PLUS_A         = 3
      integer, parameter :: PARMETIS                = 4
      integer, parameter :: MY_PERMC                = 5
      integer, parameter :: NOTRANS                 = 0 ! trans_t
      integer, parameter :: TRANS                   = 1
      integer, parameter :: CONJ                    = 2
      integer, parameter :: NOEQUIL                 = 0 ! DiagScale_t  Need?
      integer, parameter :: ROW                     = 1
      integer, parameter :: COL                     = 2
      integer, parameter :: BOTH                    = 3
      integer, parameter :: NOREFINE                = 0 ! IterRefine_t
      integer, parameter :: SINGLE                  = 1
      integer, parameter :: DOUBLE                  = 2
      integer, parameter :: EXTRA                   = 3
      integer, parameter :: LUSUP                   = 0 ! MemType  Need?
      integer, parameter :: UCOL                    = 1
      integer, parameter :: LSUB                    = 2
      integer, parameter :: USUB                    = 3
      integer, parameter :: SYSTEM                  = 0 ! LU_space_t  Need?
      integer, parameter :: USER                    = 1
      integer, parameter :: SLU_NC                  = 0 ! Stype_t
      integer, parameter :: SLU_NCP                 = 1
      integer, parameter :: SLU_NR                  = 2
      integer, parameter :: SLU_SC                  = 3
      integer, parameter :: SLU_SCP                 = 4
      integer, parameter :: SLU_SR                  = 5
      integer, parameter :: SLU_DN                  = 6
      integer, parameter :: SLU_NR_loc              = 7
      integer, parameter :: SLU_S                   = 0 ! Dtype_t
      integer, parameter :: SLU_D                   = 1
      integer, parameter :: SLU_C                   = 2
      integer, parameter :: SLU_Z                   = 3
      integer, parameter :: SLU_GE                  = 0 ! Mtype_t
      integer, parameter :: SLU_TRLU                = 1
      integer, parameter :: SLU_TRUU                = 2
      integer, parameter :: SLU_TRL                 = 3
      integer, parameter :: SLU_TRU                 = 4
      integer, parameter :: SLU_SYL                 = 5
      integer, parameter :: SLU_SYU                 = 6
      integer, parameter :: SLU_HEL                 = 7
      integer, parameter :: SLU_HEU                 = 8

      type    :: superlu_t
      integer(cptr) :: grid
      integer(cptr) :: options
      integer(cptr) :: ScalePermstruct
      integer(cptr) :: LUstruct
      integer(cptr) :: SOLVEstruct
      integer(cptr) :: A
      integer(cptr) :: stat
      end type superlu_t

      contains


      subroutine create_superlu_t( LU, nodes )
      implicit none
      integer,          intent(in) :: nodes
      type(superlu_t), intent(out) :: LU
      integer                      :: nr, nc

!     Allocate memory for superlu structures
      call f_create_gridinfo_handle( LU%grid )
      call f_create_options_handle( LU%options )
      call f_create_ScalePerm_handle( LU%ScalePermstruct )
      call f_create_LUstruct_handle( LU%LUstruct )
      call f_create_SOLVEstruct_handle( LU%SOLVEstruct )
      call f_create_SuperMatrix_handle( LU%A )
      call f_create_SuperLUStat_handle( LU%stat )

      call set_superlu_2d_grid( nodes, nr, nc )

!     Initialize the SuperLU_DIST process grid
      call f_superlu_gridinit( MPI_COMM_WORLD, nr, nc, LU%grid )

!     Set the default input options
      call f_set_default_options( LU%options )

      end subroutine create_superlu_t


      subroutine destroy_superlu_t( LU )
      implicit none
      type(superlu_t), intent(inout) :: LU
!     Release the SuperLU process grid
      call f_superlu_gridexit( LU%grid )

!     Deallocate memory of superlu structures
      call f_destroy_SuperLUStat_handle( LU%stat )
      call f_destroy_SuperMatrix_handle( LU%A )
      call f_destroy_SOLVEstruct_handle( LU%SOLVEstruct )
      call f_destroy_LUstruct_handle( LU%LUstruct )
      call f_destroy_ScalePerm_handle( LU%ScalePermstruct )
      call f_destroy_options_handle( LU%options )
      call f_destroy_gridinfo_handle( LU%grid )
      end subroutine destroy_superlu_t

      subroutine set_superlu_2d_grid( NP, NR, NC )
      implicit none
      integer,  intent(in) :: NP
      integer, intent(out) :: NR, NC
      integer              :: N1, N2

      N1 = 1
      N2 = NP

      NR = 1
      NC = NP

      do while( N1.LT.N2 )
        N1 = N1 + 1
        N2 = NP/N1
        if (MOD(NP,N1).eq.0) then
          NR = N1
          NC = N2
        else
          N2 = N2 + 1
        endif
      enddo
      end subroutine set_superlu_2d_grid

      subroutine set_superlu_options( opt, Fact, Equil, ParSymbFact,
     &                                ColPerm, RowPerm, IterRefine,
     &                                Trans, ReplaceTinyPivot,
     &                                SolveInitialized,
     &                                RefineInitialized, PrintStat )
      implicit none
      integer(cptr)     :: opt
      integer, optional :: Fact, Equil, ParSymbFact, ColPerm, RowPerm,
     &                     IterRefine, Trans, ReplaceTinyPivot,
     &                     SolveInitialized, RefineInitialized,
     &                     PrintStat
      integer           :: l_Fact, l_Equil, l_ParSymbFact, l_ColPerm,
     &                     l_RowPerm, l_IterRefine, l_Trans,
     &                     l_ReplaceTinyPivot, l_SolveInitialized,
     &                     l_RefineInitialized, l_PrintStat

      call f_get_superlu_options( opt, l_Fact, l_Equil, l_ParSymbFact,
     &                            l_ColPerm, l_RowPerm, l_IterRefine,
     &                            l_Trans, l_ReplaceTinyPivot,
     &                            l_SolveInitialized,
     &                            l_RefineInitialized,
     &                            l_PrintStat )

      if( present(Fact) )        l_Fact = Fact
      if( present(Equil) )       l_Equil = Equil
      if( present(ParSymbFact) ) l_ParSymbFact = ParSymbFact
      if( present(ColPerm) )     l_ColPerm = ColPerm
      if( present(RowPerm) )     l_RowPerm = RowPerm
      if( present(IterRefine) )  l_IterRefine = IterRefine
      if( present(Trans) )       l_Trans = Trans
      if( present(ReplaceTinyPivot) )
     &                           l_ReplaceTinyPivot = ReplaceTinyPivot
      if( present(SolveInitialized) )
     &                           l_SolveInitialized = SolveInitialized
      if( present(RefineInitialized) )
     &                           l_RefineInitialized = RefineInitialized
      if( present(PrintStat) )   l_PrintStat = PrintStat

      call f_set_superlu_options( opt, l_Fact, l_Equil, l_ParSymbFact,
     &                            l_ColPerm, l_RowPerm, l_IterRefine,
     &                            l_Trans, l_ReplaceTinyPivot,
     &                            l_SolveInitialized,
     &                            l_RefineInitialized,
     &                            l_PrintStat )
      end subroutine set_superlu_options

      subroutine get_CompRowLoc_Matrix( A, nrow, ncol, nnz_loc,
     &                                  nrow_loc, fst_row )
      implicit none
      integer(cptr)     :: A
      integer, optional :: nrow, ncol, nnz_loc, nrow_loc, fst_row
      integer           :: l_nrow, l_ncol, l_nnz_loc, l_nrow_loc,
     &                     l_fst_row

      call f_get_CompRowLoc_Matrix( A, l_nrow, l_ncol, l_nnz_loc,
     &                              l_nrow_loc, l_fst_row )

      if( present(nrow) ) nrow = l_nrow
      if( present(ncol) ) ncol = l_ncol
      if( present(nnz_loc) ) nnz_loc = l_nnz_loc
      if( present(nrow_loc) ) nrow_loc = l_nrow_loc
      if( present(fst_row) ) fst_row = l_fst_row

      end subroutine get_CompRowLoc_Matrix

      end module SUPERLU
