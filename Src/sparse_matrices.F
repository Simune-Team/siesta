      MODULE sparse_matrices
      use precision
      implicit none
!     Max. nonzero elements in a row of the hamiltonian matrix
      integer:: nh 

!     Max. number of nonzero H matrix elements    
      integer :: maxnh = 10 

      integer, pointer :: listh(:), listhptr(:),
     &                    numh(:), numhold(:)

      real(dp), pointer :: Dold(:,:), Dscf(:,:),
     &                     Eold(:,:), Escf(:,:), H(:,:), xijo(:,:),
     &                     H0(:), S(:)

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
      call de_alloc( xijo,        'xijo',        'sparseMat' )
      call de_alloc( H0,          'H0',          'sparseMat' )
      call de_alloc( S,           'S',           'sparseMat' )

      end subroutine resetSparseMatrices

      END MODULE sparse_matrices
