      MODULE sparse_matrices
      use precision
      use m_sparse_types, only: sparsity_t, sparseMatrix_t
      use m_sparse_types, only: new_sparse_matrix, release_sparse_matrix

      implicit none

!     Max. number of nonzero H matrix elements    
      integer :: maxnh = 10 

      integer, pointer :: listh(:), listhptr(:),
     &                    numh(:), numhold(:)

      real(dp), pointer :: Dold(:,:), Dscf(:,:),
     &                     Eold(:,:), Escf(:,:), H(:,:), xijo(:,:),
     &                     S(:)

      ! Pieces of H that do not depend on the SCF density matrix
      ! Formerly there was a single array H0 for this
      real(dp), pointer    :: H_kin(:) => null()

      !
      !  New interface data
      type(sparsity_t), pointer     :: sparsity => null()
      type(sparseMatrix_t)          :: H_vkb

      CONTAINS

      subroutine resetSparseMatrices( )
      use alloc,           only : de_alloc
      implicit none

      call de_alloc( listh,       'listh',       'sparseMat' )
      call de_alloc( listhptr,    'listhptr',    'sparseMat' )
      call de_alloc( numh,        'numh',        'sparseMat' )
      call de_alloc( numhold,     'numhold',     'sparseMat' )

      call de_alloc( Dold,        'Dold',        'sparseMat' )
      call de_alloc( Dscf,        'Dscf',        'sparseMat' )
      call de_alloc( Eold,        'Eold',        'sparseMat' )
      call de_alloc( Escf,        'Escf',        'sparseMat' )
      call de_alloc( H,           'H',           'sparseMat' )
      call de_alloc( H_kin,       'H_kin',       'sparseMat' )
      call release_sparse_matrix(H_vkb)
      call de_alloc( xijo,        'xijo',        'sparseMat' )
      call de_alloc( S,           'S',           'sparseMat' )

      end subroutine resetSparseMatrices

      END MODULE sparse_matrices
