MODULE sparse_matrices
  use precision
  implicit none

  ! Max. nonzero elements in a row of the hamiltonian matrix
  integer:: nh 

  ! Max. number of nonzero H matrix elements    
  integer :: maxnh = 10 



  integer, pointer :: listh(:), listhold(:), listhptr(:), listhptrold(:),  &
                      numh(:), numhold(:)

  real(dp), pointer :: Dold(:,:), Dscf(:,:), Dscfsave(:,:), Eold(:,:), &
                       Escf(:,:), H(:,:)
  real(dp), pointer :: xijo(:,:)

  real(dp), pointer :: H0(:), S(:)

END MODULE sparse_matrices
