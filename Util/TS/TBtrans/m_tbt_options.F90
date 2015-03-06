! Options for tbtrans

module m_tbt_options

  use precision, only : dp

  use m_ts_electype
  use m_ts_chem_pot

  use dictionary

  implicit none

  ! Common flags for parameters
  public
  save

  ! The temperature
  real(dp) :: kT

  ! Electrodes and different chemical potentials
  integer :: N_Elec = 0
  type(Elec), allocatable, target :: Elecs(:)
  integer :: N_mu = 0
  type(ts_mu), allocatable, target :: mus(:)

  ! Whether we should stop right after having created
  ! the Green's function files
  logical :: stop_after_GS = .false.

  ! Dictionary to contain the data saving methods
  ! Each key corresponds to some data calculation
  ! algorithm.
  ! To check whether data should be calculated do:
  ! if ( 'DOS-Gf' .in. save_DATA ) then
  !   calculate DOS of Gf
  ! end fi
  type(dict) :: save_DATA

  ! IO optimization
  ! Bigger number saves more calculations in memory
  ! But saves IO at each energy point.
!  integer :: N_io_step = 1

  ! If the energy-contour is not perfectly divisable by the number of nodes then adjust
  integer :: opt_TriMat_method = 0 ! Optimization method for determining the best tri-diagonal matrix split
  ! 0  == We optimize for speed
  ! 1  == We optimize for memory

  ! A quantity describing the accuracy of the coordinates of the 
  ! electrodes.
  ! * Should only be edited by experienced users *
  real(dp) :: Elecs_xa_EPS = 1.e-4_dp

  ! Every 5% of the calculation progress it will print an estimation
  integer :: percent_tracker = 5

#ifdef NCDF_4
  ! Save file names for data files
  character(len=400), save :: cdf_fname = ' '
  character(len=400), save :: cdf_fname_sigma = ' '
  character(len=400), save :: cdf_fname_proj = ' '
#endif

contains

  subroutine tbt_options(spin_idx, na_u, xa, lasto)

    use files, only : slabel
    use fdf
    use parallel, only : IONode
    use units, only: eV, Ang, Kelvin

    use m_tbt_save
#ifdef NCDF_4
    use m_tbt_sigma_save
#endif
    use m_tbt_diag, only : init_diag
    use m_ts_method

    use m_ts_cctype

    use m_tbt_contour

    integer, intent(in) :: spin_idx, na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: lasto(0:na_u)

    real(dp) :: Volt, rtmp
    logical :: err, ltmp
    character(len=200) :: chars
    integer :: i, j

    ! Read in the temperature
    kT = fdf_get('ElectronicTemperature',1.9e-3_dp,'Ry')

    percent_tracker = fdf_get('TBT.Progress',5)
    percent_tracker = max(1,percent_tracker)

    ! Reading the Transiesta solution method
    chars = fdf_get('TBT.SolutionMethod','BTD')
    if ( leqi(chars,'full') ) then
       !ts_method = TS_SPARSITY
       call die('Currently unsupported solution method: '//trim(chars))
    else if ( leqi(chars,'BTD').or.leqi(chars,'tri') ) then
       ts_method = TS_SPARSITY_TRI
    else
       call die('Unrecognized TBtrans solution method: '//trim(chars))
    end if

    chars = fdf_get('TS.BTD.Optimize','speed')
    chars = fdf_get('TBT.BTD.Optimize',trim(chars))
    if ( leqi(chars,'speed') ) then
       opt_TriMat_method = 0
    else if ( leqi(chars,'memory') ) then
       opt_TriMat_method = 1
    else
       call die('Could not determine flag TBT.BTD.Optimize, please &
            &see manual.')
    end if

    ! Whether we should calculate the reflectance of the junction

    ! Read in the chemical potentials
    Volt = fdf_get('TS.Voltage',0._dp,'Ry')
    Volt = fdf_get('TBT.Voltage',Volt,'Ry')
    N_mu = fdf_nmu('TBT',mus)
    if ( N_mu < 1 ) then
       N_mu = fdf_nmu('TS',mus)
    end if
    ltmp = .true.
    if ( N_mu < 1 ) then
       ltmp = .false.
       N_mu = fdffake_mu(mus,kT,Volt)
    end if
    do i = 1 , N_mu
       ! Default things that could be of importance
       if ( fdf_mu('TBT',mus(i),kT,Volt) ) then
          ! success
       else if ( fdf_mu('TS',mus(i),kT,Volt) ) then
          ! success
       else if ( ltmp ) then
          call die('Could not find chemical potential: ' &
               //trim(name(mus(i))))
       end if
    end do
#ifdef TBT_PHONON
    ! Phonon transport cannot define different chemical potentials
    ! Furthermore, they should be zero
    do i = 1 , N_mu
       if ( abs(mus(i)%mu) > 1.e-10 * eV ) then
          call die('Phonon transport does not define chemical &
               &potentials. I.e. you cannot lift the frequency spetra.')
       end if
    end do
#endif

    ! To determine the same coordinate nature of the electrodes
    Elecs_xa_EPS = fdf_get('TS.Elecs.Coord.Eps',1.e-4_dp,'Bohr')
    Elecs_xa_EPS = fdf_get('TBT.Elecs.Coord.Eps',Elecs_xa_EPS,'Bohr')

    ! detect how many electrodes we have
    N_Elec = fdf_nElec('TBT',Elecs)
    if ( N_Elec < 1 ) N_Elec = fdf_nElec('TS',Elecs)
    if ( N_Elec < 1 ) then
       ! We initialize to 2 electrodes (Left/Right)
       N_Elec = 2
       allocate(Elecs(N_Elec))
       Elecs(1)%name = 'Left'
       Elecs(1)%ID = 1
       Elecs(2)%name = 'Right'
       Elecs(2)%ID = 2
       ! if they do-not exist, the user will be told
       if ( IONode ) then
          chars = '(''tbtrans: ***'',a)'
          write(*,chars) 'No electrode names were found, &
               &default Left/Right are expected'
       end if
    end if

    ! Setup default parameters for the electrodes
    ! first electrode is the "left"
    ! last electrode is the "right"
    ! the remaining electrodes have their chemical potential at 0
    ! Currently the transport direction for all electrodes is the default
    ! We should probably warn if +2 electrodes are used and t_dir is the
    ! same for all electrodes... Then the user needs to know what (s)he is doing...
    Elecs(:)%Bulk  = fdf_get('TS.Elecs.Bulk',.true.) ! default everything to bulk electrodes
    Elecs(:)%Bulk  = fdf_get('TBT.Elecs.Bulk',Elecs(1)%Bulk)

    ! We default to not calculate the band-bottom...
    Elecs(:)%ReUseGF = fdf_get('TS.Elecs.GF.ReUse',.true.)
    Elecs(:)%ReUseGF = fdf_get('TBT.Elecs.GF.ReUse',Elecs(1)%ReUseGF)
    rtmp             = fdf_get('TS.Elecs.Eta',0.00001*eV,'Ry')
    rtmp             = fdf_get('TBT.Elecs.Eta',rtmp,'Ry')
    Elecs(:)%Eta = rtmp

    ! whether all calculations should be performed
    ! "out-of-core" i.e. whether the GF files should be created or not
    ! In tbtrans this is now defaulted to in-core
    Elecs(:)%out_of_core = fdf_get('TBT.Elecs.Out-of-core',.false.)

    do i = 1 , N_Elec
       ! Default things that could be of importance
       if ( fdf_Elec('TBT',slabel,Elecs(i),N_mu,mus) ) then
          ! success
       else if ( fdf_Elec('TS',slabel,Elecs(i),N_mu,mus,name_prefix='TBT') ) then
          ! success
       else
          call die('Could not find electrode: '//trim(name(Elecs(i))))
       end if
       ! set the placement in orbitals
       if ( Elecs(i)%idx_a < 0 ) &
            Elecs(i)%idx_a = na_u + Elecs(i)%idx_a + 1
       if ( Elecs(i)%idx_a < 1 .or. &
            na_u < Elecs(i)%idx_a ) then
          print *,Elecs(i)%idx_a,na_u
          call die("Electrode position does not exist")
       end if
       Elecs(i)%idx_o = lasto(Elecs(i)%idx_a-1)+1

       ! we need to correct the GF file name in case of
       ! single spin
       select case ( spin_idx ) 
       case ( 1 )
          ! We are using spin up
          Elecs(i)%GFfile = trim(Elecs(i)%GFfile)//'_UP'
       case ( 2 ) 
          ! We are using spin down
          Elecs(i)%GFfile = trim(Elecs(i)%GFfile)//'_DW'
       end select

       if ( (rtmp > 0._dp .and. Elecs(i)%Eta < 0._dp) .or. &
            (rtmp < 0._dp .and. Elecs(i)%Eta > 0._dp) ) then
          call die('All Eta must be either positive or negative &
               &to ensure that the retarded or advanced self-energy &
               &is exclusively calculated.')
       end if

    end do

    ! Populate the electrodes in the chemical potential type
    do i = 1 , N_Elec
       err = .true.
       do j = 1 , N_mu
          if ( associated(Elecs(i)%mu,target=mus(j)) ) then
             call chem_pot_add_Elec(mus(j),i)
             err = .false.
             exit
          end if
       end do
       if ( err ) then
          call die('We could not attribute a chemical potential &
               &to electrode: '//trim(Elecs(i)%name))
       end if
    end do

    ! check that all electrodes and chemical potentials are paired in
    ! some way.
    if ( any(mus(:)%N_El == 0) ) then
       call die('A/Some chemical potential(s) has/have not been assigned any electrodes. &
            &All chemical potentials *MUST* be assigned an electrode')
    end if

    if ( na_u <= sum(TotUsedAtoms(Elecs)) ) then
       call die('Electrodes occupy the entire structure!!!')
    end if

    ! CHECK THIS (we could allow it by only checking the difference...)
    if (  maxval(mus(:)%mu) - minval(mus(:)%mu) - abs(Volt) > 1.e-9_dp ) then
       if ( IONode ) then
          write(*,'(a)') 'Chemical potentials [eV]:'
          do i = 1 , N_Elec
             write(*,'(a,f10.5,a)') trim(Name(Elecs(i)))//' at ',Elecs(i)%mu%mu/eV,' eV'
          end do
          write(*,'(a)') 'The difference must satisfy: "max(ChemPots)-min(ChemPots) - abs(Volt) > 1e-9"'
          write(*,'(a,f10.5,a)') 'max(ChemPots) at ', maxval(mus(:)%mu)/eV,' eV'
          write(*,'(a,f10.5,a)') 'min(ChemPots) at ', minval(mus(:)%mu)/eV,' eV'
          write(*,'(a,f10.5,a)') '|V| at ', abs(Volt)/eV,' eV'
       end if
       call die('Chemical potentials are not consistent with the bias applied.')
    end if

    ! Check that the bias does not introduce a gating
    if ( any(abs(mus(:)%mu) - abs(Volt) > 1.e-9_dp ) ) then
       write(*,'(a)') 'Chemical potentials must lie in the range [-V;V] with the maximum &
            &difference being V'
       call die('Chemical potentials must not introduce consistent Ef shift to the system.')
    end if

    ! Only if the user has particular slow IO will it help 
    ! to increase this number
    ! Currently this has not been implemented, as it didn't seem to matter much
    !N_io_step = fdf_get('TBT.IO.Step',1)

    ltmp = fdf_get('TBT.R',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('T-reflect'.kv.1)
    end if

    ! Whether we should assert and calculate
    ! all transmission amplitudes
    ltmp = fdf_get('TBT.T.Elecs.All',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('T-all'.kv.1)
    end if

    ! Should we calculate DOS of Green's function
    ltmp = fdf_get('TBT.DOS.Gf', N_Elec == 1 )
    if ( ltmp ) then
       save_DATA = save_DATA // ('DOS-Gf'.kv.1)
    end if

    ! Should we calculate DOS of spectral function
    ltmp = fdf_get('TBT.DOS.A',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('DOS-A'.kv.1)
    end if

    ! Should we calculate DOS of spectral function
    ltmp = fdf_get('TBT.DOS.A.Elecs.All',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('DOS-A-all'.kv.1)
       save_DATA = save_DATA // ('DOS-A'.kv.1)
    end if

    ! Should we calculate DOS of spectral function
    ltmp = fdf_get('TBT.Current.Orb', .false. )
    if ( ltmp .and. ('DOS-A'.in.save_DATA)) then
       save_DATA = save_DATA // ('orb-current'.kv.1)
    else if ( ltmp .and. IONode ) then
       write(*,'(a,/,a)')'WARNING: Will not calculate the orbital currents, &
            &the spectral function needs to be calculated for this to &
            &apply.','Set TBT.DOS.A T to calculate orbital currents.'
    end if

    !i = fdf_get('TBT.T.Eig',0)
    if ( i > 0 ) then
       ! currently not working
       !save_DATA = save_DATA // ('T-eig'.kv.i)
    end if

    ! Will stop after creating the GF files.
    stop_after_GS = fdf_get('TBT.Elecs.GF.Only',.false.)

    call tbt_read_contour_options(N_Elec, Elecs, N_mu, mus)

    if ( IONode ) then
       write(*,7) 'Electronic temperature',kT/Kelvin,'K'
       write(*,6) 'Voltage', Volt/eV,'Volts'
       write(*,1) 'Saving DOS from Green function',('DOS-Gf'.in.save_DATA)
       if ( 'DOS-A-all' .in. save_DATA ) then
          write(*,1) 'Saving DOS from all spectral functions',.true.
       else if ( 'DOS-A' .in. save_DATA ) then
          write(*,1) 'Saving DOS from spectral functions',.true.
       else
          write(*,1) 'Saving DOS from spectral functions',.false.
       end if
       write(*,1) 'Calc. T between all electrodes',('T-all'.in.save_DATA)
       write(*,1) 'Calc. "reflection"',('T-reflect'.in.save_DATA)
       if ( spin_idx == 0 ) then
          write(*,11) 'Calculate for all spin-channels'
       else
          write(*,11) 'Only calculate for spin-channel',spin_idx
       end if
       write(*,10)'          >> Electrodes << '
       do i = 1 , size(Elecs)
          call print_settings(Elecs(i),'tbt_options')
       end do

    end if

    ! Initialize the diagonalization method.
    call init_diag( )

#ifdef NCDF_4
    call init_Sigma_options( save_DATA )
#endif
    call init_save_options( )

    if ( ('orb-current' .in.save_DATA) ) then
       ltmp = .not. fdf_get('SpinSpiral',.false.)
       ltmp = fdf_get('TBT.Symmetry.TimeReversal',ltmp)
       if ( IONode .and. ltmp ) then
          write(*,'(a,/,a)') 'WARNING: k-averaging orbital currents with &
               &time-reversal symmetry will not reproduce','the correct &
               &orbital current. Set TBT.Symmetry.TimeReversal F'
       end if
    end if


    if ( IONode ) then
       
       call print_contour_tbt_options( 'TBT' )

    end if

    ! save the used weights and energy-points.
    call io_contour_tbt(slabel)

1   format('tbt_options: ',a,t53,'=',4x,l1)
5   format('tbt_options: ',a,t53,'=',i5,a)
20  format('tbt_options: ',a,t53,'= ',i0,' -- ',i0)
6   format('tbt_options: ',a,t53,'=',f10.4,tr1,a)
7   format('tbt_options: ',a,t53,'=',f12.6,tr1,a)
8   format('tbt_options: ',a,t53,'=',f10.4)
10  format('tbt_options: ',a,t53,'=',4x,a)
11  format('tbt_options: ',a)
15  format('tbt_options: ',a,t53,'= ',i0,' x ',i0,' x ',i0)

  end subroutine tbt_options

end module m_tbt_options
