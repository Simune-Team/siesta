! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
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

      integer, public, pointer :: listh(:)=>null(), listhptr(:)=>null(),
     $                            numh(:)=>null()

      real(dp), public, pointer :: Dold(:,:)=>null(), Dscf(:,:)=>null()
      real(dp), public, pointer :: Eold(:,:)=>null(), Escf(:,:)=>null()
      real(dp), public, pointer :: Hold(:,:)=>null(), H(:,:)=>null()

      real(dp), public, pointer :: S(:)=>null()

      real(dp), public, pointer :: xijo(:,:)=>null()


      ! Pieces of H that do not depend on the SCF density matrix
      ! Formerly there was a single array H0 for this
      type(dSpData1D), public, save :: H_vkb_1D, H_kin_1D
      ! LDA+U and spin-orbit coupling Hamiltonian
      type(dSpData2D), public, save :: H_ldau_2D, H_so_2D

      !  New interface data
      type(Sparsity), public, save :: sparse_pattern
      type(dSpData1D), public, save :: S_1D
      type(dSpData2D), public, save :: DM_2D, EDM_2D, xij_2D, H_2D
      type(Fstack_Pair_Geometry_dSpData2D), public, save :: DM_history

      ! Orbital distribution descriptor
      type(OrbitalDistribution), public, save :: block_dist
      ! Always store a single distributor
      type(OrbitalDistribution), public, save :: single_dist

      public :: resetSparseMatrices

      CONTAINS

      subroutine resetSparseMatrices( )
      use alloc, only : de_alloc

      implicit none

      call delete( block_dist )
      call delete( sparse_pattern )
      nullify(numh,listhptr,listh)

      call delete( H_kin_1D )
      call delete( H_vkb_1D )
      call delete( H_ldau_2D )
      call delete( H_so_2D )

      call delete( DM_2D ) ; nullify(Dscf)
      call delete( EDM_2D ) ; nullify(Escf)
      call delete( S_1D ) ; nullify(S)
      call delete( H_2D ) ; nullify(H)
      call delete( xij_2D ) ; nullify(xijo)

      call delete( DM_history )
      
      ! Using MixH is bad as several utilities
      ! are not relying on FoX, better to leave out
      ! Note that it is nullified from start, 
      ! and hence de_alloc returns immediately
      
      call de_alloc( Hold, 'Hold', 'sparseMat' )
      call de_alloc( Dold, 'Dold', 'sparseMat' )
      call de_alloc( Eold, 'Eold', 'sparseMat' )

      end subroutine resetSparseMatrices

      END MODULE sparse_matrices
