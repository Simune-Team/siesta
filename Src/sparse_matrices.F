      MODULE sparse_matrices
      use precision
      use class_dSpData1D
      use class_dSpData2D
      use class_Sparsity
      use class_OrbitalDistribution
      use class_Fstack_Pair_Geometry_dSpData2D

      implicit none

      private

!     Max. number of nonzero H matrix elements    
      integer, public :: maxnh = 10 

      integer, public, pointer :: listh(:), listhptr(:), numh(:)

      real(dp), public, pointer :: Dold(:,:), Dscf(:,:), 
     &                             Escf(:,:), H(:,:),
     &                             xijo(:,:), S(:)

      ! Used to hold the "input" H when mixing the Hamiltonian
      real(dp), dimension(:,:), public, pointer :: Hold => null()   

      ! Pieces of H that do not depend on the SCF density matrix
      ! Formerly there was a single array H0 for this
      type(dSpData1D), public, save :: H_vkb_1D, H_kin_1D

      !  New interface data
      type(Sparsity), public, save :: sparse_pattern
      type(dSpData1D), public, save :: S_1D
      type(dSpData2D), public, save :: DM_2D, EDM_2D, xij_2D, H_2D
      type(Fstack_Pair_Geometry_dSpData2D), public, save :: DM_history

      ! Orbital distribution descriptor
      type(OrbitalDistribution), public, save :: block_dist

      public :: resetSparseMatrices

      CONTAINS

      subroutine resetSparseMatrices( )
      use alloc, only : de_alloc

      call de_alloc( Dold, 'Dold', 'state_init' )
      call de_alloc( Hold, 'Hold', 'state_init' )
      call delete( DM_2D ) ; nullify(Dscf)
      call delete( EDM_2D ) ; nullify(Escf)
      call delete( H_kin_1D )
      call delete( H_vkb_1D )
      call delete( S_1D ) ; nullify(S)
      call delete( xij_2D ) ; nullify(xijo)
      call delete( H_2D ) ; nullify(H)
      call delete( DM_history )
      call delete( sparse_pattern ) ; nullify(numh,listhptr,listh)
      call delete( block_dist )

      end subroutine resetSparseMatrices

      END MODULE sparse_matrices
