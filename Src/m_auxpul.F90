module m_auxpul

  use precision, only: dp
	implicit none

  public

  integer :: nauxpul = 1 ! Max. required size of auxiliary sparse Pulay matrices
  real(dp), allocatable   :: auxpul(:,:) ! Auxilialy matrices for Pulay mixing

end module m_auxpul





