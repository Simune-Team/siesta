MODULE m_stress
  use precision
  implicit none

  ! Constrained stress tensor
  real(dp):: cstress(3,3)

  ! Total stress tensor, includin kinetic components
  real(dp):: tstress(3,3) 

  ! Stress tensor = d_E/d_strain
  real(dp):: stress(3,3) 

  ! Local-mode part of the stress tensor
  real(dp):: stressl(3,3) 
END MODULE m_stress
