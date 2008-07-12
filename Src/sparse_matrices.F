      MODULE sparse_matrices
      use precision
      implicit none
!     Max. nonzero elements in a row of the hamiltonian matrix
      integer:: nh 

!     Max. number of nonzero H matrix elements    
      integer :: maxnh = 10 

      integer, pointer :: listh(:), listhold(:), listhptr(:),
     &                    listhptrold(:), numh(:), numhold(:)

      real(dp), pointer :: Dold(:,:), Dscf(:,:), Dscfsave(:,:),
     &                     Eold(:,:), Escf(:,:), H(:,:), xijo(:,:),
     &                     H0(:), S(:)

      CONTAINS

      subroutine resetSparseMatrices( )
      use alloc,           only : de_alloc
      implicit none
      CHARACTER(LEN=*), PARAMETER :: MYNAME =  'resetSparseMatrices'

      if (associated(listh))
     &  call de_alloc( listh, name='listh', routine=MYNAME )

      if (associated(listhold))
     &  call de_alloc( listhold, name='listhold', routine=MYNAME )

      if (associated(listhptr))
     &  call de_alloc( listhptr, name='listhptr', routine=MYNAME )

      if (associated(listhold))
     &  call de_alloc( listhold, name='listhold', routine=MYNAME )

      if (associated(listhptrold))
     &  call de_alloc( listhptrold, name='listhptrold', routine=MYNAME )

      if (associated(numh))
     &  call de_alloc( numh, name='numh', routine=MYNAME )

      if (associated(numhold))
     &  call de_alloc( numhold, name='numhold', routine=MYNAME )

      if (associated(Dold))
     &  call de_alloc( Dold, name='Dold', routine=MYNAME )

      if (associated(Dscf))
     &  call de_alloc( Dscf, name='Dscf', routine=MYNAME )

      if (associated(Dscfsave))
     &  call de_alloc( Dscfsave, name='Dscfsave', routine=MYNAME )

      if (associated(Eold))
     &  call de_alloc( Eold, name='Eold', routine=MYNAME )

      if (associated(Escf))
     &  call de_alloc( Escf, name='Escf', routine=MYNAME )

      if (associated(H))
     &  call de_alloc( H, name='H', routine=MYNAME )

      if (associated(xijo))
     &  call de_alloc( xijo, name='xijo', routine=MYNAME )

      if (associated(H0))
     &  call de_alloc( H0, name='H0', routine=MYNAME )

      if (associated(S))
     &  call de_alloc( S, name='S', routine=MYNAME )

      end subroutine resetSparseMatrices

      END MODULE sparse_matrices
