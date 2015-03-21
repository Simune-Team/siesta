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

  ! The corresponding bias for either the left electrode
  real(dp), save :: V_low = 0._dp
  real(dp), save :: V_high = 0._dp

contains

  subroutine ts_init_voltage(ucell,na_u,xa,meshG,nsm)
    use parallel, only : IONode
    use m_ts_options, only : VoltageInC, Elecs, Volt, Hartree_fname
    use units, only : eV

! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: na_u
    real(dp),      intent(in) :: xa(3,na_u)
    integer,       intent(in) :: meshG(3), nsm

    real(dp) :: tmp
    integer :: iElL, iElR

    if ( IONode ) then
       write(*,*)
       write(*,'(a,f6.3,1x,a)')'ts_voltage: Bias ', Volt/eV,'V'
    end if

    if ( ts_tdir < 1 ) then
       if ( IONode .and. len_trim(Hartree_fname) == 0 ) then
          write(*,'(a)')'ts_voltage: Lifted locally on each electrode'
       else if ( IONode .and. len_trim(Hartree_fname) > 0 ) then
          write(*,'(a)')'ts_voltage: User supplied Poisson solution'
       end if
       ! Find the lowest and highest chemical potential
       V_low = huge(1._dp)
       V_high = -huge(1._dp)
       do iElL = 1 , size(Elecs)
          V_low = min(Elecs(iElL)%mu%mu,V_low)
          V_high = max(Elecs(iElL)%mu%mu,V_high)
       end do
       if ( Volt < 0._dp ) then
          ! with a negative bias, we have to reverse
          ! the high-low
          tmp = V_high
          V_high = V_low
          V_low = tmp
       end if

       return

    end if

    ! set the left chemical potential
    call get_elec_indices(na_u, xa, iElL, iElR)
    V_low = Elecs(iElL)%mu%mu
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
    use m_ts_options, only : Hartree_fname
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
    if ( ts_tdir > 0 ) then
       call ts_ramp_elec(ucell,ntpl,Vscf)
#ifdef NCDF
    else if ( len_trim(Hartree_fname) > 0 ) then
       call ts_ncdf_Voltage(Hartree_fname,'V',ntpl,Vscf)
#endif
    else
       call ts_elec_only(ntpl,Vscf)
    end if

    call timer('ts_volt',2)

  end subroutine ts_voltage


  subroutine ts_ramp_elec(ucell, ntpl, Vscf)
    use precision,    only : grid_p
    use m_ts_options, only : Volt
    use m_mesh_node,  only : meshl, offset_i
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
    real(dp) :: dF

    ! field in [0;end]: v = e*x = f*index
    dF = Volt / real(right_elec_mesh_idx - left_elec_mesh_idx,dp)

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
                Vscf(imesh) = Vscf(imesh) + V_low - dF*idT
             end do
          end do
       end do

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
                Vscf(imesh) = Vscf(imesh) + V_low - dF*idT
             end do
          end do
       end do
       
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
                Vscf(imesh) = Vscf(imesh) + V_low - dF*idT
             end do
          end do
       end do

    end if

  end subroutine ts_ramp_elec

  subroutine ts_elec_only(ntpl, Vscf)
    use precision,    only : grid_p
    use m_ts_options, only : N_Elec, Elecs
    use m_mesh_node,  only : meshl, offset_r, dMesh, dL
    use m_geom_box, only : voxel_in_box_delta
! ***********************
! * INPUT variables     *
! ***********************
    integer,       intent(in) :: ntpl
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(ntpl)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: ix, iy, iz, iEl, imesh
    real(dp) :: llZ(3), llYZ(3), ll(3)
#ifdef TRANSIESTA_BOX
    integer, allocatable :: n_V(:)

    allocate(n_V(N_Elec))
    n_V = 0
#endif

    ! We do a loop in the local grid
    imesh = 0
    do iz = 0 , meshl(3) - 1
       llZ(:) = offset_r(:) + iz*dL(:,3)
       do iy = 0 , meshl(2) - 1
          llYZ(:) = iy*dL(:,2) + llZ(:)
          do ix = 0 , meshl(1) - 1
             ! The lower-left corner of the current box
             ll = ix*dL(:,1) + llYZ
             imesh = imesh + 1

             ! Count entries in each geometric object
             do iEl = 1 , N_Elec
                if ( voxel_in_box_delta(Elecs(iEl)%box, ll, dMesh) ) then
#ifdef TRANSIESTA_BOX
                   n_V(iEl) = n_V(iEl) + 1
#endif
                   Vscf(imesh) = Vscf(imesh) + Elecs(iEl)%mu%mu
                end if
             end do

          end do
       end do
    end do

#ifdef TRANSIESTA_BOX
    print '(10(tr1,i8))',n_V,ntpl
#endif

    if ( imesh /= ntpl ) then
       call die('Electrode bias did not loop on all points. Something &
            &needs to be checked.')
    end if

  end subroutine ts_elec_only

  subroutine get_elec_indices(na_u, xa, iElL, iElR)
    use m_ts_electype
    use m_ts_options, only : Elecs
    
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
    do i = Elecs(1)%idx_a , Elecs(1)%idx_a + TotUsedAtoms(Elecs(1)) - 1
       ! Correct for rotated unitcells TODO
       if ( abs(xa(ts_tdir,i)) < tmp ) then
          tmp = abs(xa(ts_tdir,i))
       end if
    end do
    do i = Elecs(2)%idx_a , Elecs(2)%idx_a + TotUsedAtoms(Elecs(2)) - 1
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
    do i = Elecs(iElL)%idx_a , &
         Elecs(iElL)%idx_a + TotUsedAtoms(Elecs(iElL)) - 1
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
    do i = Elecs(iElR)%idx_a , &
         Elecs(iElR)%idx_a + TotUsedAtoms(Elecs(iElR)) - 1
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
       write(*,'(a,f6.3,1x,a)')'ts_voltage: Bias @bottom ', V_low/eV,'V'
       write(*,'(a,3(f6.3,a))')'ts_voltage: In unit cell direction = {', &
            vcdir(1),',',vcdir(2),',',vcdir(3),'}'
    end if

  end subroutine print_ts_voltage

#ifdef NCDF
  ! Read in a potential file from a NetCDF file.
  ! Thus the potential landscape can be fully customized by the user.
  ! We note that the potential landscape need only be calculated
  ! for one V, direct interpolation is possible as 
  ! the solution to the Poisson equation is linearly dependent on the BC
  subroutine ts_ncdf_voltage(fname,V_name, npt,V)
    use precision, only: grid_p
    use nf_ncdf
    use dictionary
    use variable
    
    character(len=*), intent(in) :: fname, V_name
    integer, intent(in) :: npt
    real(grid_p), intent(inout) :: V(npt)
    
    type(hNCDF) :: ncdf
    type(dict) :: atts
    real(grid_p) :: Vmin, Vmax, fact
    real(grid_p), allocatable :: tmpV(:)

    ! All nodes should allocate an auxilliary grid
    allocate(tmpV(npt))

    ! Open the file
    call ncdf_open(ncdf,fname)

    ! Read in the grid
    call cdf_r_grid(ncdf,trim(V_name),lnpt,tmpV)

    ! retrieve the max min from the file
    call ncdf_get_var(ncdf,trim(V_name)//'min',Vmin)
    call ncdf_get_var(ncdf,trim(V_name)//'max',Vmax)

    call ncdf_close(ncdf)

    ! Correct the limits so that we align to the current potential
    fact = ( V_high - V_low ) / ( Vmax - Vmin ) 

!$OMP parallel workshare default(shared), firstprivate(fact)
    V = V + tmpV * fact
!$OMP end parallel workshare
    
    deallocate(tmpV)

  end subroutine ts_ncdf_voltage

#endif

end module m_ts_voltage
