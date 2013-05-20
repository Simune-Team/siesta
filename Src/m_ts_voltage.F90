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
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
module m_ts_voltage

! Module for containing various ways of implementing the voltage in a
! TranSIESTA run.
! As a standard method TranSIESTA applies the voltage ramp across the
! entire unit cell.
! As a future progress the voltage ramp might be situated in the contact region
! only. This makes more physical sense as the voltage drop does not occur
! in the electrode regions.
!
! Created and copyrighted by: Nick Papior Andersen, 2012
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.
!
  
  use precision, only : dp
  
  implicit none
  
  private

  ! The idea is to have sub routines in this module to do
  ! various voltage layouts
  public :: ts_voltage
  public :: ts_init_voltage

  ! The voltage direction (currently this module relies on the voltage direction
  ! to be only in the z-direction)
  integer, parameter :: V_DIR = 3

  ! we set the left/right indices for the input of the voltage ramp.
  ! This will put the voltage in between the electrodes instead of 
  ! in the entire cell.
  ! It should leverage some iterations as the voltage drop ensures a 
  ! faster density convergence
  integer, save :: left_elec_mesh_idx = 0
  integer, save :: right_elec_mesh_idx = huge(1)

contains

  subroutine ts_init_voltage(cell,na_u,xa,meshG,nsm)
    use m_ts_options, only : VoltageInC
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: na_u
    real(dp),      intent(in) :: xa(3,na_u)
    integer,       intent(in) :: meshG(3), nsm

    call print_ts_voltage(cell)

    if ( VoltageInC ) then
       ! Find the electrode mesh sets
       call init_elec_indices(cell, meshG, nsm, na_u, xa)
    end if

  end subroutine ts_init_voltage

  subroutine ts_voltage(cell,meshG,nsm,Vscf)
    use precision,    only : grid_p
    use m_ts_options, only : VoltageInC
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(*)
    call timer('ts_volt',1)

    ! Here one should implement different calls...

    if ( VoltageInC ) then
       ! Voltage drop in between the electrodes
       call ts_ramp_elec(cell,meshG,nsm,Vscf)
    else
       ! Voltage drop in the entire cell
       call ts_ramp_cell(cell,meshG,nsm,Vscf)
    end if

    call timer('ts_volt',2)
  end subroutine ts_voltage

  subroutine ts_ramp_cell(cell, meshG, nsm, Vscf)
    use precision,    only : grid_p
    use parallel,     only : IONode
    use mesh,         only : meshLim
    use m_ts_options, only : VoltFDF, VoltL
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(*)

! ***********************
! * LOCAL variables     *
! ***********************
    integer                   :: ix, i1, i2, i3, i10, i20, i30, imesh
    integer                   :: meshl(3), Xoffset, Yoffset, Zoffset
    real(dp)                  :: Lvc, dLvc, dF
    real(dp)                  :: dot
    external                  :: dot

    Lvc = sqrt(dot(cell(1,V_DIR),cell(1,V_DIR),3))

    dLvc = Lvc/max( meshG(V_DIR), 1 ) !

! field in [0;Lvc]: v = e*x = f*index
    dF = -VoltFDF*dLvc/Lvc

! Find quantities in mesh coordinates
    meshl(1) = (meshLim(2,1) - meshLim(1,1)+1)*nsm
    meshl(2) = (meshLim(2,2) - meshLim(1,2)+1)*nsm
    meshl(3) = (meshLim(2,3) - meshLim(1,3)+1)*nsm

! Calculate starting point for grid
    Xoffset = (meshLim(1,1)-1)*nsm
    Yoffset = (meshLim(1,2)-1)*nsm
    Zoffset = (meshLim(1,3)-1)*nsm

! Add the electric field potential to the input potential
    imesh = 0
    i30 = Zoffset - 1
    do i3 = 0,meshl(3)-1
       i30 = i30 + 1
       i20 = Yoffset - 1
       do i2 = 0,meshl(2)-1
          i20 = i20 + 1
          i10 = Xoffset - 1
          do i1 = 0,meshl(1)-1
             i10   = i10 + 1
             imesh = imesh + 1
             Vscf(imesh) = Vscf(imesh) + VoltL + dF*i30
          enddo
       enddo
    enddo

  end subroutine ts_ramp_cell

  subroutine ts_ramp_elec(cell, meshG, nsm, Vscf)
    use precision,    only : grid_p
    use parallel,     only : IONode
    use mesh,         only : meshLim
    use m_ts_options, only : VoltFDF, VoltL, VoltR
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(*)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i1, i2, i3, idz, imesh
    integer  :: meshl(3), Zoffset
    real(dp) :: Lvc, dLvc, dF
    real(dp) :: dot
    external :: dot

    Lvc = sqrt(dot(cell(1,V_DIR),cell(1,V_DIR),3))

    dLvc = Lvc/max( meshG(V_DIR), 1 )

    ! For the electrode the distance is only in between the indices
    Lvc = dLvc * (right_elec_mesh_idx - left_elec_mesh_idx)

    ! field in [0;Lvc]: v = e*x = f*index
    dF = VoltFDF*dLvc/Lvc

    ! Find quantities in mesh coordinates
    meshl(1) = (meshLim(2,1) - meshLim(1,1)+1)*nsm
    meshl(2) = (meshLim(2,2) - meshLim(1,2)+1)*nsm
    meshl(3) = (meshLim(2,3) - meshLim(1,3)+1)*nsm

    ! Calculate starting point for grid
    Zoffset = (meshLim(1,3)-1)*nsm

    ! Add the electric field potential to the input potential
    imesh = 0
    do i3 = Zoffset+1,Zoffset+meshl(3)

       ! Check the region
       if ( i3 <= left_elec_mesh_idx ) then
          ! We are on the left hand side, hence
          ! we have the left \mu
          idz = 0
       else if ( right_elec_mesh_idx < i3 ) then
          ! We are on the right hand side, hence
          ! we have the right \mu
          idz = right_elec_mesh_idx - left_elec_mesh_idx
       else
          ! We are in between the electrodes
          idz = i3 - left_elec_mesh_idx
       end if

       do i2 = 0,meshl(2)-1
          do i1 = 0,meshl(1)-1
             imesh = imesh + 1
             Vscf(imesh) = Vscf(imesh) + VoltL - dF*idz
          enddo
       enddo
    enddo

  end subroutine ts_ramp_elec

  subroutine init_elec_indices(cell, meshG, nsm, na_u, xa)
    use parallel,     only : IONode
    use units,        only : Ang
    use m_ts_options, only : na_BufL, na_BufR
    use m_ts_electype
    use m_ts_options, only : ElLeft, ElRight
    
! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: cell(3,3)
    integer,  intent(in) :: meshG(3), nsm, na_u
    real(dp), intent(in) :: xa(3,na_u)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i, iz
    real(dp) :: Lvc, dLvc
    real(dp) :: left_z_max, right_z_min
    real(dp) :: ddleft, ddright

    ! We require that the cell is deterministic in the z-direction
    Lvc = abs(cell(V_DIR,V_DIR))

    left_z_max  = -huge(1._dp)
    right_z_min = huge(1._dp)
    do i = na_BufL + 1 , na_BufL + TotUsedAtoms(ElLeft)
       if ( left_z_max < xa(V_DIR,i) ) then
          left_z_max = xa(V_DIR,i) + 0.25_dp ! We add 0.25 Bohr for a small distance to the electrode
       end if
    end do
    do i = na_u - na_BufR - TotUsedAtoms(ElRight) + 1 , na_u - na_BufR
       if ( right_z_min > xa(V_DIR,i) ) then
          right_z_min = xa(V_DIR,i) - 0.25_dp ! We add 0.25 Bohr for a small distance to the electrode
       end if
    end do

    ddleft  = huge(1._dp)
    ddright = huge(1._dp)

    ! Initialize the indices if it does not exist
    left_elec_mesh_idx  = 1
    right_elec_mesh_idx = meshG(V_DIR)
    
    ! The distance step in the z-direction
    dLvc = Lvc/max( meshG(V_DIR), 1 ) !
    do iz = 0 , meshG(V_DIR) - 1
       if ( abs(dLvc * iz - left_z_max) < ddleft ) then
          ddleft = abs(dLvc * iz - left_z_max)
          left_elec_mesh_idx = iz + 1
       end if
       if ( abs(dLvc * iz - right_z_min) < ddright ) then
          ddright = abs(dLvc * iz - right_z_min)
          right_elec_mesh_idx = iz + 1
       end if
    end do

    ! We correct for the case of users having put the electrode in 
    ! in-correct order in the cell....
    iz = 0
    if ( left_elec_mesh_idx > right_elec_mesh_idx ) then
       if ( IONode ) then
          write(*,'(1x,a)') 'WARNING: The voltage ramp could not be placed in &
               &between the electrodes, be sure to have the atomic &
               &coordinates in ascending order and starting from &
               & (0.,0.,0.).'
       end if
       left_elec_mesh_idx = 1
       right_elec_mesh_idx = meshG(V_DIR)
    end if
    if ( left_elec_mesh_idx >= meshG(V_DIR) ) then
       if ( IONode ) then
          write(*,'(1x,a)') 'WARNING: The voltage ramp could not be placed in &
               &between the electrodes, be sure to have the atomic &
               &coordinates in ascending order and starting from &
               & (0.,0.,0.).'
       end if
       left_elec_mesh_idx = 1
       right_elec_mesh_idx = meshG(V_DIR)
    end if

    if ( left_elec_mesh_idx /= 1 .and. &
         right_elec_mesh_idx /= meshG(V_DIR) .and. IONode ) then
       write(*,'(a,2(f9.3,a))') 'ts_voltage: Ramp placed between the &
            &cell z-coordinates (Ang): ',(left_elec_mesh_idx-1)*dLvc/Ang, &
            ' and ',right_elec_mesh_idx*dLvc/Ang
    end if
    
  end subroutine init_elec_indices

! Print out the voltage direction dependent on the cell parameters.

  subroutine print_ts_voltage(cell)
    use parallel,     only : IONode
    use units,        only : eV
    use m_ts_options, only : VoltFDF
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: cell(3,3)

! ***********************
! * LOCAL variables     *
! ***********************
    integer                   :: i
    real(dp)                  :: Lvc, vcdir(3)

    real(dp)                  :: dot
    external                  :: dot
    
    Lvc = sqrt(dot(cell(1,V_DIR),cell(1,V_DIR),3))
    do i=1,3
       vcdir(i) = cell(i,V_DIR)/Lvc
    end do

    if(IONode) then
       write(*,*)
       write(*,'(a,f6.3,1x,a)')'ts_voltage: Bias ', VoltFDF/eV,'V'
       write(*,'(a,3(f6.3,a))')'ts_voltage: In unit cell direction = {', &
            vcdir(1),',',vcdir(2),',',vcdir(3),'}'
    end if
  end subroutine print_ts_voltage

end module m_ts_voltage
