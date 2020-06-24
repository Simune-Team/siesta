!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

! Created by Nick Papior in June 2020.

!> Module intented for handling dipoles in the system
!!
!! This module does 2 things:
!! 1. Returns the center of a dipole which defaults to the
!!    center of the system in case the user does not know where
!!    it is.
!! 2. Calculate the dipole based on a given center.
module dipole_m

  implicit none

  private
  public :: get_dipole_center
  public :: calc_dipole

contains

  !< Return the dipole center
  !!
  !! If the user does not specify the dipole center
  !! the dipole center will be the center of the system (not center-of-mass).
  !!
  !! However, in some skewed molecules where one knows the dipole is
  !! not located at the center of the system one should specify the
  !! center of the dipole through this block: Dipole.Center
  subroutine get_dipole_center(ucell, na_u, xa, x0)
    use precision, only: dp
    use sys, only: die
    use fdf

    integer, intent(in) :: na_u
    real(dp), intent(in) :: ucell(3,3), xa(3,na_u)
    real(dp), intent(out) :: x0(3)

    ! Local variables
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=132) :: length_unit
    real(dp) :: cfactor
    integer :: ia

    x0(:) = 0.0_dp

    if ( fdf_block("Dipole.Center", bfdf) ) then

      ! TODO consider allowing atomic indices
      ! as input. Then the center could be calculated
      ! with respect to these coordinates. Here
      ! multiple values of the same atom *could* be allowed.
      ! Say a single atom above a hexagon?

      if ( .not. fdf_bline(bfdf,pline) ) then
        call die("Dipole.Center could not find center")
      end if

      if ( .not. fdf_bmatch(pline,"vvvn") ) then
        call die("Dipole.Center both a coordinate and unit *must* be present")
      end if

      length_unit = fdf_bnames(pline, 1)
      cfactor = fdf_convfac(length_unit, "Bohr")

      ! Retrieve value
      ! This will not work for MD simulations (unless
      ! the dipole center is fixed!)
      do ia = 1,3
        x0(ia) = fdf_bvalues(pline,ia) * cfactor
      end do

      call fdf_bclose(bfdf)

    else

      do ia = 1, na_u
        x0(1:3) = x0(1:3) + xa(1:3,ia)
      end do
      ! Average
      do ia = 1, 3
        x0(ia) = x0(ia) / na_u
      end do

    end if

  end subroutine get_dipole_center

  !< Calculate dipole for a given center of the dipole
  subroutine calc_dipole( cell, dvol, ntm, ntml, nsm, drho, X0, dipole)
! ********************************************************************
! Finds the electric dipole
! Written by J.M.Soler. July 1997.
! Modified for distributed drho matrix using a 2D grid of processors.
! Routine now is based on intrinsic structure of grid distribution
! in order to calculate position of grid points in space for local
! matrix. J.D.Gale March 1999.
! *********** INPUT **************************************************
! real*8  cell(3,3)     : Unit cell vectors
! real*8  dvol          : Voxel volume
! integer ntm(3)        : Global number of divisions of each lattice vector
! integer ntml(3)       : Local number of divisions of each lattice vector
! integer nsm           : Number of sub-points for each mesh point
! real*?  drho(*ntml)   : Minus neutral charge density at mesh points
! real*8  X0(3)         : Origin in cartesian coordinates of the dipole
! *********** OUTPUT *************************************************
! real*8 dipole(3)   : Electric dipole
! *********** UNITS **************************************************
! cell   in atomic units (Bohr)
! dvol   in atomic units (Bohr**3)
! drho   in atomic units (electrons/Bohr**3)
! X0     in atomic units (Bohr)
! dipole in atomic units (electrons*Bohr)
! ********************************************************************
!
!  Modules
!
    use precision,  only : dp, grid_p
    use sys,        only : die
    use mesh,       only : meshLim
#ifdef MPI
    use mpi_siesta
#endif

    integer, intent(in) :: ntm(3), ntml(3), nsm
    real(grid_p), intent(in) :: drho(ntml(1),ntml(2),ntml(3))
    real(dp), intent(in) :: cell(3,3), dvol, X0(3)
    real(dp), intent(out) :: dipole(3)

    external :: reclat

    ! Internal variables and arrays
    integer :: I, I1, I2, I3, I10, I20, I30
    integer :: MG1, MG2, MG3, Xoffset, Yoffset, Zoffset
    real(dp) :: D(3), DD(3), DX(3,2:3), Rcell(3,3), X0L(3)
#ifdef MPI
    integer :: MPIerror
#endif

    ! Assign local variables
    MG1 = ntm(1)
    MG2 = ntm(2)
    MG3 = ntm(3)

    do I = 1 , 3
      DD(I) = 1._dp / ntm(I)
    end do

    ! Find reciprocal cell vectors (without the factor 2*pi)
    call reclat(cell, Rcell, 0 )

    ! Find origin in lattice coordinates
    do I = 1,3
      X0L(I) = X0(1)*Rcell(1,I) + X0(2)*Rcell(2,I) + X0(3)*Rcell(3,I)
    end do

    ! Initialize dipole
    dipole(1:3) = 0.0_dp

    ! Calculate starting point for grid
    Xoffset = (meshLim(1,1)-1)*nsm
    Yoffset = (meshLim(1,2)-1)*nsm
    Zoffset = (meshLim(1,3)-1)*nsm

    ! Find dipole by direct integration allowing for block distributed
    ! structure of drho
    I30 = Zoffset - 1
    do I3 = 1 , ntml(3)
      I30 = I30 + 1
      D(3) = DD(3) * I30 - X0L(3)
      if ( D(3) < - 0.5_dp ) then
        D(3) = D(3) + 1.0_dp
      else if ( D(3) > + 0.5_dp ) then
        D(3) = D(3) - 1.0_dp
      end if
      DX(1,3) = cell(1,3)*D(3)
      DX(2,3) = cell(2,3)*D(3)
      DX(3,3) = cell(3,3)*D(3)
      I20 = Yoffset - 1
      do I2 = 1 , ntml(2)
        I20 = I20 + 1
        D(2) = DD(2) * I20 - X0L(2)
        if ( D(2) < - 0.5_dp ) then
          D(2) = D(2) + 1.0_dp
        else if ( D(2) > + 0.5_dp ) then
          D(2) = D(2) - 1.0_dp
        end if
        DX(1,2) = cell(1,2)*D(2) + DX(1,3)
        DX(2,2) = cell(2,2)*D(2) + DX(2,3)
        DX(3,2) = cell(3,2)*D(2) + DX(3,3)
        I10 = Xoffset - 1
        do I1 = 1 , ntml(1)
          I10 = I10 + 1
          D(1) = DD(1) * I10 - X0L(1)
          if ( D(1) < - 0.5_dp ) then
            D(1) = D(1) + 1.0_dp
          else if ( D(1) > + 0.5_dp ) then
            D(1) = D(1) - 1.0_dp
          end if
          dipole(1) = dipole(1) - (cell(1,1)*D(1) + DX(1,2)) * drho(I1,I2,I3)
          dipole(2) = dipole(2) - (cell(2,1)*D(1) + DX(2,2)) * drho(I1,I2,I3)
          dipole(3) = dipole(3) - (cell(3,1)*D(1) + DX(3,2)) * drho(I1,I2,I3)
        end do
      end do
    end do

    ! Correct for volume
    dipole(1) = dipole(1) * dvol
    dipole(2) = dipole(2) * dvol
    dipole(3) = dipole(3) * dvol

#ifdef MPI
    call MPI_AllReduce(dipole, d, 3, MPI_double_precision, MPI_sum, &
        MPI_Comm_World, MPIerror)
    dipole(1:3) = d(1:3)
#endif

  end subroutine calc_dipole

end module dipole_m
