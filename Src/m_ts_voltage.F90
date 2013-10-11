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
  use m_ts_tdir

  implicit none
  
  private

  ! The idea is to have sub routines in this module to do
  ! various voltage layouts
  public :: ts_voltage
  public :: ts_init_voltage
  public :: ts_VH_fix

  ! we set the left/right indices for the input of the voltage ramp.
  ! This will put the voltage in between the electrodes instead of 
  ! in the entire cell.
  ! It should leverage some iterations as the voltage drop ensures a 
  ! faster density convergence
  integer, save :: left_elec_mesh_idx = 0
  integer, save :: right_elec_mesh_idx = huge(1)

contains

  subroutine ts_init_voltage(ucell,na_u,xa,meshG,nsm)
    use m_ts_options, only : VoltageInC
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: na_u
    real(dp),      intent(in) :: xa(3,na_u)
    integer,       intent(in) :: meshG(3), nsm

    call print_ts_voltage(ucell)

    if ( VoltageInC ) then
       ! Find the electrode mesh sets
       call init_elec_indices(ucell, meshG, nsm, na_u, xa)
    end if

  end subroutine ts_init_voltage

  subroutine ts_voltage(ucell,meshG,nsm,Vscf)
    use precision,    only : grid_p
    use m_ts_options, only : VoltageInC
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(*)
    call timer('ts_volt',1)

    ! Here one should implement different calls...

    if ( VoltageInC ) then
       ! Voltage drop in between the electrodes
       call ts_ramp_elec(ucell,meshG,nsm,Vscf)
    else
       ! Voltage drop in the entire cell
       call ts_ramp_cell(ucell,meshG,nsm,Vscf)
    end if

    call timer('ts_volt',2)
  end subroutine ts_voltage

  subroutine ts_ramp_cell(ucell, meshG, nsm, Vscf)
    use precision,    only : grid_p
    use parallel,     only : IONode
    use mesh,         only : meshLim
    use m_ts_options, only : VoltFDF, VoltL
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(*)

! ***********************
! * LOCAL variables     *
! ***********************
    integer                   :: i1, i2, i3, imesh
    integer, target           :: i10, i20, i30
    integer, pointer          :: iT
    integer                   :: meshl(3), Xoffset, Yoffset, Zoffset
    real(dp)                  :: Lvc, dLvc, dF
    real(dp)                  :: dot
    external                  :: dot

    Lvc = sqrt(dot(ucell(1,ts_tdir),ucell(1,ts_tdir),3))

    dLvc = Lvc/max( meshG(ts_tdir), 1 ) !

    ! field in [0;Lvc]: v = e*x = f*index
    dF = -VoltFDF*dLvc/Lvc

    ! Set up counter
    if ( ts_tdir == 1 ) then
       iT => i10
    else if ( ts_tdir == 2 ) then
       iT => i20
    else
       iT => i30
    end if
 
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
             Vscf(imesh) = Vscf(imesh) + VoltL + dF*iT
          enddo
       enddo
    enddo

  end subroutine ts_ramp_cell

  subroutine ts_ramp_elec(ucell, meshG, nsm, Vscf)
    use precision,    only : grid_p
    use parallel,     only : IONode
    use mesh,         only : meshLim
    use m_ts_options, only : VoltFDF, VoltL, VoltR
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)
    integer,       intent(in) :: meshG(3), nsm
! ***********************
! * OUTPUT variables    *
! ***********************
    real(grid_p), intent(inout) :: Vscf(*)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i1, i2, i3, idT, imesh
    integer  :: meshl(3), Xoffset, Yoffset, Zoffset
    real(dp) :: Lvc, dLvc, dF
    real(dp) :: dot
    external :: dot

    Lvc = sqrt(dot(ucell(1,ts_tdir),ucell(1,ts_tdir),3))

    dLvc = Lvc/max( meshG(ts_tdir), 1 )

    ! For the electrode the distance is only in between the indices
    Lvc = dLvc * (right_elec_mesh_idx - left_elec_mesh_idx)

    ! field in [0;Lvc]: v = e*x = f*index
    dF = VoltFDF*dLvc/Lvc

    ! Find quantities in mesh coordinates
    meshl(1) = (meshLim(2,1) - meshLim(1,1)+1)*nsm
    meshl(2) = (meshLim(2,2) - meshLim(1,2)+1)*nsm
    meshl(3) = (meshLim(2,3) - meshLim(1,3)+1)*nsm

    ! Calculate starting point for grid
    Xoffset = (meshLim(1,1)-1)*nsm
    Yoffset = (meshLim(1,2)-1)*nsm
    Zoffset = (meshLim(1,3)-1)*nsm

    if ( ts_tdir == 1 ) then

       ! Add the electric field potential to the input potential
       imesh = 0
       do i3 = 1,meshl(3)
          do i2 = 1,meshl(2)
             do i1 = Xoffset+1,Xoffset+meshl(1)
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
                Vscf(imesh) = Vscf(imesh) + VoltL - dF*idT
             enddo
          enddo
       enddo

    else if ( ts_tdir == 2 ) then

       ! Add the electric field potential to the input potential
       imesh = 0
       do i3 = 1,meshl(3)
          do i2 = Yoffset+1,Yoffset+meshl(2)
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
                Vscf(imesh) = Vscf(imesh) + VoltL - dF*idT
             enddo
          enddo
       enddo
       
    else
       ! Add the electric field potential to the input potential
       imesh = 0
       do i3 = Zoffset+1,Zoffset+meshl(3)

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
                Vscf(imesh) = Vscf(imesh) + VoltL - dF*idT
             enddo
          enddo
       enddo

    end if

  end subroutine ts_ramp_elec

  subroutine init_elec_indices(ucell, meshG, nsm, na_u, xa)
    use parallel,     only : IONode
    use units,        only : Ang
    use m_ts_options, only : na_BufL, na_BufR
    use m_ts_electype
    use m_ts_options, only : Elecs
    
! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: ucell(3,3)
    integer,  intent(in) :: meshG(3), nsm, na_u
    real(dp), intent(in) :: xa(3,na_u)

! ***********************
! * LOCAL variables     *
! ***********************
    integer  :: i, it
    real(dp) :: Lvc, dLvc
    real(dp) :: left_t_max, right_t_min
    real(dp) :: ddleft, ddright
    real(dp) :: ElecL(3), ElecR(3)
    real(dp) :: dot
    external :: dot

    if ( size(Elecs) > 2 ) call die('Not fully implemented')
    Lvc = sqrt(dot(ucell(1,ts_tdir),ucell(1,ts_tdir),3))

    left_t_max  = -huge(1._dp)
    right_t_min = huge(1._dp)
    do i = na_BufL + 1 , na_BufL + TotUsedAtoms(Elecs(1))
       if ( left_t_max < xa(ts_tdir,i) ) then
          left_t_max = xa(ts_tdir,i) + 0.25_dp ! We add 0.25 Bohr for a small distance to the electrode
       end if
    end do
    do i = na_u - na_BufR - TotUsedAtoms(Elecs(2)) + 1 , na_u - na_BufR
       if ( right_t_min > xa(ts_tdir,i) ) then
          right_t_min = xa(ts_tdir,i) - 0.25_dp ! We add 0.25 Bohr for a small distance to the electrode
       end if
    end do

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

  subroutine print_ts_voltage(ucell)
    use parallel,     only : IONode
    use units,        only : eV
    use m_ts_options, only : VoltFDF
! ***********************
! * INPUT variables     *
! ***********************
    real(dp),      intent(in) :: ucell(3,3)

! ***********************
! * LOCAL variables     *
! ***********************
    integer                   :: i
    real(dp)                  :: Lvc, vcdir(3)

    real(dp)                  :: dot
    external                  :: dot
    
    Lvc = sqrt(dot(ucell(1,ts_tdir),ucell(1,ts_tdir),3))
    do i=1,3
       vcdir(i) = ucell(i,ts_tdir)/Lvc
    end do

    if(IONode) then
       write(*,*)
       write(*,'(a,f6.3,1x,a)')'ts_voltage: Bias ', VoltFDF/eV,'V'
       write(*,'(a,3(f6.3,a))')'ts_voltage: In unit cell direction = {', &
            vcdir(1),',',vcdir(2),',',vcdir(3),'}'
    end if
  end subroutine print_ts_voltage



! Fix the voltage in the mesh
  subroutine ts_VH_fix(  mesh, nsm, v )
!
!  Modules
!
    use precision, only : dp, grid_p
    use sys, only : die
    use parallel, only : ProcessorY, Node, Nodes
#ifdef MPI
    use mpi_siesta, only : MPI_AllReduce, MPI_Sum
    use mpi_siesta, only : MPI_Integer, MPI_Comm_World
    use mpi_siesta, only : DAT_double => MPI_double_precision
    use parallelsubs, only : HowManyMeshPerNode
#endif
    
    implicit          none
    real(grid_p)      v(*)
    integer           mesh(3), nsm

! Internal variables
    integer           i1, i2, i3, imesh, ntemp
    integer           nlp, meshl(3)
    integer           ProcessorZ, BlockSizeY, BlockSizeZ, Yoffset
    integer           Zoffset, Py, Pz, meshnsm(3)
    integer, target  :: i10, i20, i30
    integer, pointer :: iT
    integer           NRemY, NRemZ
#ifdef MPI
    integer           MPIerror, npl
#endif
    real(dp)          vav, vtot, temp

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE TSVHfix' )
#endif

    call timer('TSVHFix',1)

    ! Set up counter
    if ( ts_tdir == 1 ) then
       iT => i10
    else if ( ts_tdir == 2 ) then
       iT => i20
    else
       iT => i30
    end if

! Find local number of mesh points
    meshnsm(1) = mesh(1)/nsm
    meshnsm(2) = mesh(2)/nsm
    meshnsm(3) = mesh(3)/nsm
#ifdef MPI
!! print *, "N:",Node, "TSVHFix-- mesh:", mesh
    call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
    meshl(1) = meshl(1)*nsm
    meshl(2) = meshl(2)*nsm
    meshl(3) = meshl(3)*nsm
#else
    meshl(1) = mesh(1)
    meshl(2) = mesh(2)
    meshl(3) = mesh(3)
#endif

! Check that ProcessorY is a factor of the number of processors
    if (mod(Nodes,ProcessorY).gt.0) then
       call die('ERROR: ProcessorY must be a factor of the &
            &number of processors!')
    endif
    ProcessorZ = Nodes/ProcessorY

! Calculate blocking sizes
    BlockSizeY = (meshnsm(2)/ProcessorY)*nsm
    NRemY = (mesh(2) - ProcessorY*BlockSizeY)/nsm
    BlockSizeZ = (meshnsm(3)/ProcessorZ)*nsm
    NRemZ = (mesh(3) - ProcessorZ*BlockSizeZ)/nsm

! Calculate coordinates of current node in processor grid
    Py = (Node/ProcessorZ)+1
    Pz = Node - (Py - 1)*ProcessorZ + 1

! Calculate starting point for grid
    Yoffset = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY)
    Zoffset = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ)

    vtot = 0._dp
    nlp  = meshl(1)*meshl(2)*meshl(3)

    ! Test whether we should do anything (note that iT => [i10|i20|i30]):
    i10 = 0
    i20 = Yoffset - 1
    i30 = Zoffset - 1
    if ( iT <= 0 ) then
       imesh = 0
       i30 = Zoffset - 1
       do i3 = 0,meshl(3)-1
          i30 = i30 + 1
          i20 = Yoffset - 1
          do i2 = 0,meshl(2)-1
             i20 = i20 + 1
             do i10 = 0,meshl(1)-1
                imesh = imesh + 1
                if (iT.eq.0) then
                   vtot = vtot + v(imesh)
                endif
             enddo
          enddo
       enddo
    end if

#ifdef MPI
    temp=0.d0
    call MPI_AllReduce(vtot,temp,1,DAT_double,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    vtot = temp
    ntemp=0
    call MPI_AllReduce(nlp,ntemp,1,MPI_integer,MPI_Sum, &
         MPI_Comm_World,MPIerror)
    nlp = ntemp

#endif
    vav = vtot/real(nlp,dp)

    imesh = 0
    do i30 = 1 , meshl(3)
       do i20 = 1 , meshl(2)
          do i10 = 1 , meshl(1)
             imesh = imesh + 1
             v(imesh) = v(imesh) - vav
          enddo
       enddo
    enddo

    call timer('TSVHFix',2)

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS TSVHfix' )
#endif

  end subroutine ts_VH_fix


end module m_ts_voltage
