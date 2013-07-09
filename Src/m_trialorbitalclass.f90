module trialorbitalclass

use precision,       only:dp

implicit none

type trialorbital
  real(dp),dimension(3) :: center ! Projection function center in 
                                  !   cristallographic coordinates relative 
                                  !   to the direct lattice vectors.
  real(dp),dimension(3) :: zaxis  ! Defines the axis from which the polar
                                  !   angle theta in spherical polar coordinates
                                  !   is measured.
                                  !   Default: (0.0 0.0 1.0)
  real(dp),dimension(3) :: xaxis  ! Defines the axis from which the azimuthal 
                                  !   angle phi in spherical coordinates is 
                                  !   measured.
                                  !   Must be orthogonal to z-axis.
  real(dp),dimension(3) :: yaxis  ! Angular momentum y-axis
  real(dp)              :: zovera ! z/a, diffusivity, spread. 
                                  !   Read from the nnkp file in Ang^-1
                                  !   Transformed later to Ang^-1
  integer               :: r      ! Radial quantum number
  integer               :: l      ! Angular momentum
  integer               :: mr     ! z-projection quantum number
  real(dp)              :: rcut   ! Siesta's cut-off radius: Bohr
  integer               :: lmax   ! Maximum total angular momentum
end type

! Cut-off radii in units of \alpha^-1
real(dp),parameter,dimension(3) ::                           &
                           cutoffs = (/6.934_dp,18.87_dp,35.44_dp/)
! Squared norm tolerance governing the cut-off radii
real(dp),parameter :: T = 0.0001_dp

endmodule trialorbitalclass
