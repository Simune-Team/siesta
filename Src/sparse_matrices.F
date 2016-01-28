! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE sparse_matrices
      use precision
      implicit none

!     Max. number of nonzero H matrix elements    
      integer :: maxnh = 10 

      integer, pointer :: listh(:), listhptr(:),
     &                    numh(:), numhold(:)

      real(dp), pointer :: Dold(:,:), Dscf(:,:),
     &                     Eold(:,:), Escf(:,:), H(:,:), xijo(:,:),
     &                     S(:)

      ! Used to hold the "input" H when mixing the Hamiltonian
      real(dp), dimension(:,:), pointer :: Hold => null()   

      ! Pieces of H that do not depend on the SCF density matrix
      ! Formerly there was a single array H0 for this
      real(dp), pointer :: H_kin(:) => null()
      real(dp), pointer :: H_vkb(:) => null()

      CONTAINS

      subroutine resetSparseMatrices( )
      use alloc,           only : de_alloc
!      use siesta_options, only  : MixH
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
      ! Remove pernicious dependency on siesta_options via MixH
      if (associated(Hold)) then
         call de_alloc( Hold,        'Hold',        'sparseMat' )
      endif
      call de_alloc( H_kin,       'H_kin',       'sparseMat' )
      call de_alloc( H_vkb,       'H_vkb',       'sparseMat' )
      call de_alloc( xijo,        'xijo',        'sparseMat' )
      call de_alloc( S,           'S',           'sparseMat' )

      end subroutine resetSparseMatrices

      END MODULE sparse_matrices
