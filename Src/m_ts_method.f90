module m_ts_method
  
  implicit none

  public
  save
  
  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA.
  integer, parameter :: TS_SPARSITY = 1

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA as well 
  ! as the heavily optimized tri-diagonalization
  integer, parameter :: TS_SPARSITY_TRI = 2


  ! The default solution method (it will be reset
  ! after option reading)
  integer :: ts_method = TS_SPARSITY_TRI

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

contains

! Function for calculating which region in the scattering matrix the
! hamiltonian elements adhere.
  pure function get_scat_region(io,noL,jo, noR,no_u)
    integer, intent(in) :: io, noL, jo, noR, no_u
    integer :: get_scat_region

    if ( io < 1 .and. jo < 1 ) then
       get_scat_region = 1 ! Left buffer
    else if ( (io < 1 .and. jo <= noL) .or. &
         (io <= noL .and. jo < 1) ) then
       get_scat_region = 2 ! Left buffer / left electrode
    else if ( (1 <= io .and. 1 <= jo) .and. &
         (io <= noL .and. jo <= noL) ) then
       get_scat_region = 3 ! Left electrode
    else if ( ((1 <= io .and. noL < jo) .and. &
         (io <= noL .and. jo <= no_u-noR)) .or. &
         ((noL < io .and. 1 <= jo) .and. &
         (io <= no_u-noR .and. jo <= noL)) ) then
       get_scat_region = 4 ! Left electrode / contact region
    else if ( (noL < io .and. noL < jo) .and. &
         (io <= no_u-noR .and. jo <= no_u-noR) ) then
       get_scat_region = 5 ! Contact region
    else if ( ((noL < io .and. no_u-noR < jo) .and. &
         (io <= no_u-noR .and. jo <= no_u)) .or. &
         ((no_u-noR < io .and. noL < jo) .and. &
         (io <= no_u .and. jo <= no_u-noR)) ) then
       get_scat_region = 6 ! Contact region / right electrode
    else if ( (no_u-noR < io .and. no_u-noR < jo) .and. &
         (io <= no_u .and. jo <= no_u) ) then
       get_scat_region = 7 ! Right electrode
    else if ( ((no_u-noR < io .and. no_u < jo) .and. &
         io <= no_u) .or. &
         ((no_u < io .and. no_u-noR < jo) .and. &
         jo <= no_u) ) then
       get_scat_region = 8 ! Right electrode / right buffer
    else if ( no_u < io .and. no_u < jo ) then
       get_scat_region = 9 ! Right buffer
    else
       get_scat_region = 0 ! Everything else (LE/RE ... etc.)
    end if

  end function get_scat_region

end module m_ts_method
  
