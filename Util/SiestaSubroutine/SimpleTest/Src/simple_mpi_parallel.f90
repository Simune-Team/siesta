! ---
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program simple

! A very simple driver for Siesta-as-subroutine (or siesta-as-server)
! This version uses MPI and siesta as a subroutine. It must be compiled
! together with siesta.

  use mpi
  use fsiesta

  implicit none
  integer,parameter:: dp = kind(1.d0)

  integer,parameter :: na = 3
  integer :: error, ia, myNode=0
  real(dp):: e, fa(3,na), xa(3,na)

  data xa / 0.0, 0.0, 0.0, &
            0.7, 0.7, 0.0, &
           -0.7, 0.7, 0.0 /

! Initialize MPI and get my node's index
  call MPI_Init( error )
  call MPI_Comm_Rank( MPI_Comm_World, myNode, error )

! Set physical units
  call siesta_units( 'Ang', 'eV' )

! Launch a siesta process using all available MPI processes
  call siesta_launch( 'h2o' )
  if (myNode==0) print*, 'siesta launched'

! Find forces
  call siesta_forces( 'h2o', na, xa, energy=e, fa=fa )
  if (myNode==0) &
    print'(/,a,/,(3f12.6,3x,3f12.6))', 'xa, fa =', (xa(:,ia),fa(:,ia),ia=1,na)

! Find forces for another geometry
  xa(1,1) = 0.1
  call siesta_forces( 'h2o', na, xa, energy=e, fa=fa )
  if (myNode==0) &
    print'(/,a,/,(3f12.6,3x,3f12.6))', 'xa, fa =', (xa(:,ia),fa(:,ia),ia=1,na)

! Quit siesta process
  call siesta_quit( 'h2o' )

! Finalize MPI
  call MPI_Finalize( error )

end program simple

