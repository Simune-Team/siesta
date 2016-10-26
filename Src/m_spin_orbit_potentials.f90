module m_spin_orbit_potentials
  public :: valid_spin_orbit_potentials

CONTAINS

  ! Checks whether a PSML file contains the proper semilocal
  ! potentials to construct the spin-orbit part
  
  function valid_spin_orbit_potentials(ps) result(p)
    use m_psml
    
    type(ps_t), intent(in) :: ps
    logical                :: p

    integer :: nso, nlj
    
    nso = ps_Number_Of_Potentials(ps,SET_SO)
    nlj = ps_Number_Of_Potentials(ps,SET_LJ)

    ! Crude test for now
    
    p = (nso > 0) .or. (nlj > 0)
    
  end function valid_spin_orbit_potentials

end module m_spin_orbit_potentials
    
