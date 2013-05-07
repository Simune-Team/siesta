module m_ts_method
  
  implicit none

  public
  save
  
  ! This is the original transiesta
  integer, parameter :: TS_ORIGINAL = 0

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA.
  integer, parameter :: TS_SPARSITY = 1

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA as well 
  ! as the heavily optimized tri-diagonalization
  integer, parameter :: TS_SPARSITY_TRI = 2


  ! The default solution method (it will be reset
  ! after option reading)
  integer :: ts_method = TS_ORIGINAL

  ! Logical variable that describes the solution method on
  ! LEFT-RIGHT-EQUILIBRIUM contour points.
  ! This is realized by the fact that for:
  !    UseBulk .and. UpdateDMCR
  ! the needed part of the GF is only the C...C regions:
  !
  !  -------------------------------------
  !  | L...L | L...C   0     0   |   0   |
  !  | C...L | C...C C...C C...C |   0   |
  !  |   0   | C...C C...C C...C |   0   |
  !  |   0   | C...C C...C C...C | C...R |
  !  |   0   |   0     0   R...C | R...R |
  !  -------------------------------------
  ! 
  ! This means we can solve the following instead:
  ! G_F^{-1} G_F I_P = I \times I_P,
  ! where I_P:
  !  ---------------------
  !  |   0     0     0   |
  !  |   1     0     0   |
  !  |   0     1     0   |
  !  |   0     0     1   |
  !  |   0     0     0   |
  !  ---------------------
  ! Note, that this can ONLY be used in EQUI contour points.
  ! In principle we can obtain the EXACT size of the problem
  ! For very large electrodes. This could come in handy.
  logical, save :: GF_INV_EQUI_PART = .false.

end module m_ts_method
  
