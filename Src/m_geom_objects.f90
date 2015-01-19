!
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.
!

! This module collects all the geometrical objects in the sub modules
!    m_geom_plane
!    m_geom_square
!    m_geom_coord

module m_geom_objects

  use m_geom_aux
  use m_geom_plane
  use m_geom_square
  use m_geom_coord
  use m_geom_box
  
  implicit none

  private :: dp

  interface voxel_in
     module procedure voxel_in_plane_delta
     module procedure voxel_in_plane_gauss
     module procedure voxel_in_plane_exp
     module procedure voxel_in_square_delta
     module procedure voxel_in_square_gauss
     module procedure voxel_in_square_exp
     module procedure voxel_in_box_delta
     module procedure voxel_in_coord_exp
  end interface 

  interface voxel_val
     module procedure voxel_val_plane_delta
     module procedure voxel_val_plane_gauss
     module procedure voxel_val_plane_exp
     module procedure voxel_val_square_delta
     module procedure voxel_val_square_gauss
     module procedure voxel_val_square_exp
     module procedure voxel_val_box_delta
     module procedure voxel_val_coord_exp
  end interface 

  interface fgeo_read
     module procedure fgeo_read_plane_delta
     module procedure fgeo_read_plane_gauss
     module procedure fgeo_read_plane_exp
     module procedure fgeo_read_square_delta
     module procedure fgeo_read_square_gauss
     module procedure fgeo_read_square_exp
     module procedure fgeo_read_box_delta
     module procedure fgeo_read_coord_exp
  end interface 

end module m_geom_objects
