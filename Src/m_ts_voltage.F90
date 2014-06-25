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
  
  use precision, only : dp
  use m_ts_tdir

  implicit none
  
  private

  ! The idea is to have sub routines in this module to do
  ! various voltage layouts
  public :: ts_voltage
  public :: ts_init_voltage

  ! we set the left/right indices for the input of the voltage ramp.
  ! This will put the voltage in between the electrodes instead of 
  ! in the entire cell.
  ! It should leverage some iterations as the voltage drop ensures a 
  ! faster density convergence
  integer, save :: left_elec_mesh_idx = 0
  integer, save :: right_elec_mesh_idx = huge(1)

  ! The corresponding bias' for the two different
  ! electrodes
  real(dp), save :: left_V = 0._dp

contains

  subroutine ts_init_voltage(ucell,na_u,xa,meshG,nsm)
    use m_ts_options, only : VoltageInC, Elecs

! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: na_u
    real(dp),      intent(in) :: xa(3,na_u)
    integer,       intent(in) :: meshG(3), nsm

    integer :: iElL, iElR

    if ( ts_tdir < 1 ) then
       call die('Non-determined bias direction not implemented')
    end if

    ! set the left chemical potential
    call get_elec_indices(na_u, xa, iElL, iElR)
    left_V = Elecs(iElL)%mu%mu
    call print_ts_voltage(ucell)

    if ( VoltageInC ) then
       ! Find the electrode mesh sets
       call init_elec_indices(ucell, meshG, nsm, na_u, xa)
    else
       ! Simulate the electrodes at the ends
       ! This leverages a double routine
       left_elec_mesh_idx  = 1
       right_elec_mesh_idx = meshG(ts_tdir)
    end if

  end subroutine ts_init_voltage

  subroutine ts_voltage(ucell, ntpl, Vscf)
    use precision,    only : grid_p
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: ntpl
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(ntpl)

    call timer('ts_volt',1)

    ! Voltage drop in between the electrodes
    ! The indices for the full cell is set
    ! correctly to not have two routines doing the
    ! same
    call ts_ramp_elec(ucell,ntpl,Vscf)

    call timer('ts_volt',2)

  end subroutine ts_voltage


  subroutine ts_ramp_elec(ucell, ntpl, Vscf)
    use intrinsic_missing, only : VNORM
    use precision,    only : grid_p
    use m_ts_options, only : Volt
    use m_ts_mesh,    only : meshl, offset_i, dMesh
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: ntpl
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(ntpl)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i1, i2, i3, idT, imesh
    real(dp) :: Lvc, dLvc, dF

    dLvc = dMesh(ts_tdir)
    Lvc = VNORM(ucell(:,ts_tdir))

    ! For the electrode the distance is only in between the indices
    Lvc = dLvc * (right_elec_mesh_idx - left_elec_mesh_idx)

    ! field in [0;Lvc]: v = e*x = f*index
    dF = Volt * dLvc / Lvc

    ! Find quantities in mesh coordinates
    if ( meshl(1) * meshl(2) * meshl(3) /= ntpl ) &
         call die('ERROR: Vscf size not correct')

    ! Add the electric field potential to the input potential
    imesh = 0
    if ( ts_tdir == 1 ) then

       do i3 = 1,meshl(3)
          do i2 = 1,meshl(2)
             do i1 = offset_i(1)+1,offset_i(1)+meshl(1)
                ! Check the region
                if ( i1 <= left_elec_mesh_idx ) then
                   ! We are on the left hand side, hence
                   ! we have the left \mu
                   idT = 0
                else if ( right_elec_mesh_idx < i1 ) then
                   ! We are on the right hand side, hence
                   ! we have the right \mu
                   idT = right_elec_mesh_idx - left_elec_mesh_idx
                else
                   ! We are in between the electrodes
                   idT = i1 - left_elec_mesh_idx
                end if
                
                imesh = imesh + 1
                Vscf(imesh) = Vscf(imesh) + left_V - dF*idT
             enddo
          enddo
       enddo

    else if ( ts_tdir == 2 ) then

       do i3 = 1,meshl(3)
          do i2 = offset_i(2)+1,offset_i(2)+meshl(2)
             ! Check the region
             if ( i2 <= left_elec_mesh_idx ) then
                ! We are on the left hand side, hence
                ! we have the left \mu
                idT = 0
             else if ( right_elec_mesh_idx < i2 ) then
                ! We are on the right hand side, hence
                ! we have the right \mu
                idT = right_elec_mesh_idx - left_elec_mesh_idx
             else
                ! We are in between the electrodes
                idT = i2 - left_elec_mesh_idx
             end if
             do i1 = 1,meshl(1)
                imesh = imesh + 1
                Vscf(imesh) = Vscf(imesh) + left_V - dF*idT
             enddo
          enddo
       enddo
       
    else

       do i3 = offset_i(3)+1,offset_i(3)+meshl(3)

          ! Check the region
          if ( i3 <= left_elec_mesh_idx ) then
             ! We are on the left hand side, hence
             ! we have the left \mu
             idT = 0
          else if ( right_elec_mesh_idx < i3 ) then
             ! We are on the right hand side, hence
             ! we have the right \mu
             idT = right_elec_mesh_idx - left_elec_mesh_idx
          else
             ! We are in between the electrodes
             idT = i3 - left_elec_mesh_idx
          end if

          do i2 = 1,meshl(2)
             do i1 = 1,meshl(1)
                imesh = imesh + 1
                Vscf(imesh) = Vscf(imesh) + left_V - dF*idT
             enddo
          enddo
       enddo

    end if

  end subroutine ts_ramp_elec

  subroutine get_elec_indices(na_u, xa, iElL, iElR)
    use m_ts_electype
    use m_ts_options, only : N_Elec, Elecs
    
! ***********************
! * INPUT variables     *
! ***********************
    integer,  intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)

! ***********************
! * OUTPUT variables    *
! ***********************
    integer, intent(out) :: iElL, iElR

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i
    real(dp) :: tmp

    ! find "lower" electrode
    tmp = huge(1._dp)
    iElL = 1
    iElR = 2
    do i = Elecs(1)%idx_na , Elecs(1)%idx_na + TotUsedAtoms(Elecs(1)) - 1
       ! Correct for rotated unitcells TODO
       if ( abs(xa(ts_tdir,i)) < tmp ) then
          tmp = abs(xa(ts_tdir,i))
       end if
    end do
    do i = Elecs(2)%idx_na , Elecs(2)%idx_na + TotUsedAtoms(Elecs(2)) - 1
       ! Correct for rotated unitcells TODO
       if ( abs(xa(ts_tdir,i)) < tmp ) then
          tmp = abs(xa(ts_tdir,i))
          iElL = 2
          iElR = 1
       end if
    end do

  end subroutine get_elec_indices

  subroutine init_elec_indices(ucell, meshG, nsm, na_u, xa)
    use intrinsic_missing, only : VNORM
    use parallel,     only : IONode
    use units,        only : Ang
    use m_ts_electype
    use m_ts_options, only : N_Elec, Elecs
    
! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: ucell(3,3)
    integer,  intent(in) :: meshG(3), nsm, na_u
    real(dp), intent(in) :: xa(3,na_u)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i, it, iElL, iElR
    real(dp) :: Lvc, dLvc
    real(dp) :: left_t_max, right_t_min
    real(dp) :: ddleft, ddright
    real(dp) :: ElecL(3), ElecR(3)

    if ( N_Elec > 2 ) call die('Not fully implemented, only non-bias with N-electrode')
    Lvc = VNORM(ucell(:,ts_tdir))

    ! get the left/right electrodes
    call get_elec_indices(na_u,xa,iElL,iElR)

    if ( Elecs(iElL)%Bulk ) then
       left_t_max  = -huge(1._dp)
    else
       left_t_max  =  huge(1._dp)
    end if
    if ( Elecs(iElR)%Bulk ) then
       right_t_min =  huge(1._dp)
    else
       right_t_min = -huge(1._dp)
    end if
    do i = Elecs(iElL)%idx_na , &
         Elecs(iElL)%idx_na + TotUsedAtoms(Elecs(iElL)) - 1
       if ( Elecs(iElL)%Bulk ) then
          ! Correct for rotated unitcells TODO
          if ( left_t_max < xa(ts_tdir,i) ) then
             left_t_max = xa(ts_tdir,i)
          end if
       else
          if ( xa(ts_tdir,i) < left_t_max ) then
             left_t_max = xa(ts_tdir,i)
          end if
       end if
    end do
    left_t_max = left_t_max + 0.25_dp ! We add 0.25 Bohr for a small distance to the electrode
    do i = Elecs(iElR)%idx_na , &
         Elecs(iElR)%idx_na + TotUsedAtoms(Elecs(iElR)) - 1
       if ( Elecs(iElL)%Bulk ) then
          ! Correct for rotated unitcells TODO
          if ( xa(ts_tdir,i) < right_t_min ) then
             right_t_min = xa(ts_tdir,i)
          end if
       else
          if ( right_t_min < xa(ts_tdir,i) ) then
             right_t_min = xa(ts_tdir,i)
          end if
       end if
    end do
    right_t_min = right_t_min - 0.25_dp ! We add 0.25 Bohr for a small distance to the electrode

    ddleft  = huge(1._dp)
    ddright = huge(1._dp)

    ! Initialize the indices if it does not exist
    left_elec_mesh_idx  = 1
    right_elec_mesh_idx = meshG(ts_tdir)
    
    ! The distance step in the t-direction
    dLvc = Lvc/max( meshG(ts_tdir), 1 ) !
    do it = 0 , meshG(ts_tdir) - 1
       if ( abs(dLvc * it - left_t_max) < ddleft ) then
          ddleft = abs(dLvc * it - left_t_max)
          left_elec_mesh_idx = it + 1
       end if
       if ( abs(dLvc * it - right_t_min) < ddright ) then
          ddright = abs(dLvc * it - right_t_min)
          right_elec_mesh_idx = it + 1
       end if
    end do

    ! We correct for the case of users having put the electrode in 
    ! in-correct order in the cell....
    it = 0
    if ( left_elec_mesh_idx > right_elec_mesh_idx ) then
       if ( IONode ) then
          write(*,'(1x,a)') 'WARNING: The voltage ramp could not be placed in &
               &between the electrodes, be sure to have the atomic &
               &coordinates in ascending order and starting from &
               & (0.,0.,0.).'
       end if
       left_elec_mesh_idx  = 1
       right_elec_mesh_idx = meshG(ts_tdir)
    end if
    if ( left_elec_mesh_idx >= meshG(ts_tdir) ) then
       if ( IONode ) then
          write(*,'(1x,a)') 'WARNING: The voltage ramp could not be placed in &
               &between the electrodes, be sure to have the atomic &
               &coordinates in ascending order and starting from &
               & (0.,0.,0.).'
       end if
       left_elec_mesh_idx  = 1
       right_elec_mesh_idx = meshG(ts_tdir)
    end if

    ElecL = 0._dp
    ElecR = 0._dp
    ElecL(ts_tdir) = (left_elec_mesh_idx-1)*dLvc/Ang
    ElecR(ts_tdir) = right_elec_mesh_idx*dLvc/Ang

    if ( left_elec_mesh_idx /= 1 .and. &
         right_elec_mesh_idx /= meshG(ts_tdir) .and. IONode ) then
       write(*,'(a,/,a,3(f9.3,tr1),a,3(tr1,f9.3),a)') 'ts_voltage: Ramp placed between the &
            &cell coordinates (Ang):',' {',ElecL,'} to {',ElecR,' }'
    end if
    
  end subroutine init_elec_indices

! Print out the voltage direction dependent on the cell parameters.

  subroutine print_ts_voltage( ucell )
    use intrinsic_missing, only : VNORM
    use parallel,     only : IONode
    use units,        only : eV
    use m_ts_options, only : Volt
! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: ucell(3,3)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i
    real(dp) :: Lvc, vcdir(3)

    Lvc = VNORM(ucell(:,ts_tdir))
    do i = 1 , 3
       vcdir(i) = ucell(i,ts_tdir)/Lvc
    end do

    if ( IONode ) then
       write(*,*)
       write(*,'(a,f6.3,1x,a)')'ts_voltage: Bias ', Volt/eV,'V'
       write(*,'(a,f6.3,1x,a)')'ts_voltage: Bias @bottom ', left_V/eV,'V'
       write(*,'(a,3(f6.3,a))')'ts_voltage: In unit cell direction = {', &
            vcdir(1),',',vcdir(2),',',vcdir(3),'}'
    end if

  end subroutine print_ts_voltage

end module m_ts_voltage
