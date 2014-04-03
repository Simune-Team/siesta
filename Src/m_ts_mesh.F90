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
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
module m_ts_mesh

! Module for retaining information about the mesh.
! It is used by the Hartree module and the bias module.
!
! Created and copyrighted by: Nick Papior Andersen, 2014
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.

  use precision, only : dp

  implicit none

  ! The offset for the current node
  real(dp), save :: offset_r(3) = 0._dp
  ! the voxel-vectors for each sub-mesh element
  real(dp), save :: dL(3,3) = 0._dp
  ! the voxel length along each direction
  real(dp), save :: dMesh(3) = 0._dp

  ! offsets for the local node
  integer, save :: meshl(3) = 0
  integer, save :: offset_i(3) = 0

contains

  subroutine ts_init_mesh(ucell,meshG,meshLim,nsm)
    use intrinsic_missing, only : VNORM
    use parallel, only : Node, Nodes, IONode, ProcessorY

    ! The unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of mesh divisions of each lattice vector
    integer, intent(in) :: meshG(3), meshLim(2,3)
    ! Number of fine points per big point (see iogrid_netcdf)
    integer, intent(in) :: nsm 

    ! Number of big division points
    integer :: nm(3)
    ! Processor specifics
    integer :: ProcessorZ, blocY, blocZ, nremY, nremZ
    ! dimension tracking of the divisions
    integer :: iniX, iniY, iniZ, dimX, dimY, dimZ
    ! Local node dimensionality of the grid
    integer :: ldimX, ldimY, ldimZ
    ! Loop stuff
    integer :: node_Y, node_Z, cur_Node
    integer :: i
    
    ! We calculate the spacing in each direction
    ! Notice that we now have:
    ! dL(1,1) = dX for stepping in the x-direction
    ! dL(1,2) = dX for stepping in the y-direction
    ! dL(1,3) = dX for stepping in the z-direction
    do i = 1 , 3 
       ldimX = max(meshG(i),1)
       ! The dimension stepping in each direction.
       dL(:,i) = ucell(:,i) / ldimX
       ! The mesh box-size
       dMesh(i) = VNORM(ucell(:,i)) / ldimX
    end do

    ! For nodes == 1 we have no offset
    ! (also some of the arrays are not initialized, which
    !  could lead to errors)
    if ( Nodes == 1 ) then
       meshl = meshG ! same as meshLim(2,:) * nsm
       return
    end if
       
    ! Now we need to calculate the offset of the local node

    ! Calculate the number of big-points
    nm(1:3) = meshG(1:3) / nsm

    ! In order to check that we do it correctly...
    dimX  = nm(1) * nsm
    ldimX = dimX

    ! Calculate number of z-processor divisions
    ProcessorZ = Nodes / ProcessorY
    
    ! Calculate the block-division in the Y-direction
    blocY = nm(2) / ProcessorY
    ! Calculate the block-division in the Z-direction
    blocZ = nm(3) / ProcessorZ
    ! If there are any remaining mesh-points, they will be added one 
    ! to each of the processors, in sequence
    nremY = nm(2) - blocY*ProcessorY
    nremZ = nm(3) - blocZ*ProcessorZ

    ! Initialize the loop construct variables
    cur_Node = 0
    ! Notice that we start from zero as we would like to calculate the
    ! offset
    iniX = 0
    iniY = 0
    do node_Y = 1, ProcessorY

       ! Initialize the dimension size in the Y-direction
       dimY = blocY
       ! The if there were too many points to be perfectly divisable
       ! we will add them.
       if ( node_Y <= nremY ) dimY = dimY + 1
       ! Extend into the mesh size
       dimY = dimY * nsm

       ! Initialize the Z-division (notice we index from zero
       ! as we would like to calculate the offset)
       iniZ = 0
       do node_Z = 1, ProcessorZ
          dimZ = blocZ
          if ( node_Z <= nremZ ) dimZ = dimZ + 1
          dimZ = dimZ * nsm ! For fine points
           
          if ( Node == cur_Node ) then
             ! Calculate the offset in the [YZ]-direction for the processor
             ! We know that iniX == 0, so no need, but we have it for
             ! consistency
             offset_r(:) = iniX * dL(:,1) + iniY * dL(:,2) + iniZ * dL(:,3)
             ldimY = dimY
             ldimZ = dimZ
          end if
          
          iniZ = iniZ + dimZ
          cur_Node = cur_Node + 1
       end do
       iniY = iniY + dimY
    end do

    ! Find quantities in mesh coordinates
    meshl(1) = (meshLim(2,1) - meshLim(1,1)+1)*nsm
    meshl(2) = (meshLim(2,2) - meshLim(1,2)+1)*nsm
    meshl(3) = (meshLim(2,3) - meshLim(1,3)+1)*nsm

    ! Calculate starting point for grid
    offset_i(1) = (meshLim(1,1)-1)*nsm
    offset_i(2) = (meshLim(1,2)-1)*nsm
    offset_i(3) = (meshLim(1,3)-1)*nsm

  end subroutine ts_init_mesh

end module m_ts_mesh

