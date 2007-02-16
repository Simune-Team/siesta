MODULE sparse_matrices
  use precision
  implicit none
  integer, pointer :: listh(:), listhold(:), listhptr(:), listhptrold(:),  &
                      numh(:), numhold(:)

  real(dp), pointer :: Dold(:,:), Dscf(:,:), Dscfsave(:,:), Eold(:,:), &
                       Escf(:,:), H(:,:)

  real(dp), pointer :: H0(:), S(:)

END MODULE sparse_matrices
