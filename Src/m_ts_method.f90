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
  
