      MODULE sparse_matrices
      use precision
      use class_dSpArr2D
      use class_Sparsity
      use class_OrbitalDistribution
      use class_Fstack_Pair_Geometry_dSpArr2D

      implicit none


!     Max. number of nonzero H matrix elements    
      integer :: maxnh = 10 

      integer, pointer :: listh(:), listhptr(:),
     &                    numh(:), numhold(:)

      real(dp), pointer :: Dold(:,:), Dscf(:,:), 
     &                     Escf(:,:), H(:,:),
     &                     xijo(:,:), S(:)

      ! Used to hold the "input" H when mixing the Hamiltonian
      real(dp), dimension(:,:), pointer :: Hold => null()   

      ! Pieces of H that do not depend on the SCF density matrix
      ! Formerly there was a single array H0 for this
      real(dp), pointer    :: H_kin(:) => null()

      !
      !  New interface data
      type(Sparsity), save                :: sparse_pattern
      type(dSpArr2D), save                :: H_vkb, DM, EDM
      type(Fstack_Pair_Geometry_dSpArr2D), save :: DM_history

      ! Orbital distribution descriptor
      type(OrbitalDistribution), save     :: block_dist

      CONTAINS

      subroutine resetSparseMatrices( )
      use alloc,           only : de_alloc
      implicit none

      call de_alloc( listh,       'listh',       'sparseMat' )
      call de_alloc( listhptr,    'listhptr',    'sparseMat' )
      call de_alloc( numh,        'numh',        'sparseMat' )
      call de_alloc( numhold,     'numhold',     'sparseMat' )

      call de_alloc( Dold,        'Dold',        'sparseMat' )
!      call de_alloc( Dscf,        'Dscf',        'sparseMat' )
      call delete(DM)
!      call de_alloc( Escf,        'Escf',        'sparseMat' )
      call delete(EDM)
      call de_alloc( H,           'H',           'sparseMat' )
      call de_alloc( H_kin,       'H_kin',       'sparseMat' )
      call delete(H_vkb)
      call de_alloc( xijo,        'xijo',        'sparseMat' )
      call de_alloc( S,           'S',           'sparseMat' )
      call delete(sparse_pattern)
      call delete(DM_history)
      call delete(block_dist)

      end subroutine resetSparseMatrices

      END MODULE sparse_matrices
