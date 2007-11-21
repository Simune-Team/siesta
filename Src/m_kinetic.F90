module m_kinetic

  use precision, only: dp
	implicit none

  public

  real(dp):: vn         ! Velocity (time derivative) of the Nose thermostat
  real(dp):: vpr        ! Velocity (time derivative) of the PR variables

  real(dp):: tempion=0.0_dp ! Ionic temperature

  real(dp):: kn         ! Kinetic energy of the Nose' thermostat
  real(dp):: kpr        ! Kinetic energy of the Parrinello-Rahman variables

end module m_kinetic



