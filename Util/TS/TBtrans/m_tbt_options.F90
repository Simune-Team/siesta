! Options for tbtrans

module m_tbt_options

  use precision, only : dp

  use m_ts_tdir, only: ts_tdir, ts_tidx

  use m_ts_electype
  use m_ts_chem_pot

  use dictionary

  implicit none

  ! Common flags for parameters
  public
  save

  ! The standard name_prefix
#ifdef TBT_PHONON
  character(len=*), parameter :: name_prefix = 'PHT'
#else
  character(len=*), parameter :: name_prefix = 'TBT'
#endif

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

  ! Number of eigenchannels to calculate
  integer :: N_eigen = 0

  ! If the energy-contour is not perfectly divisable by the number of nodes then adjust
  integer :: BTD_method = 0 ! Optimization method for determining the best tri-diagonal matrix split
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
  character(len=250) :: cdf_fname = ' '
  character(len=250) :: cdf_fname_sigma = ' '
  character(len=250) :: cdf_fname_proj = ' '
#endif


  ! List of private formats for printing information
  character(len=*), parameter, private :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'
  character(len=*), parameter, private :: f10='(''tbt: '',a,t53,''='',tr4,a)'
  character(len=*), parameter, private :: f11='(''tbt: '',a)'
  character(len=*), parameter, private :: f12='(''tbt: '',a,t53,''='',tr2,i0)'
  character(len=*), parameter, private :: f5 ='(''tbt: '',a,t53,''='',i5,a)'
  character(len=*), parameter, private :: f20='(''tbt: '',a,t53,''='',i0,'' -- '',i0)'
  character(len=*), parameter, private :: f6 ='(''tbt: '',a,t53,''='',f10.4,tr1,a)'
  character(len=*), parameter, private :: f7 ='(''tbt: '',a,t53,''='',f12.6,tr1,a)'
  character(len=*), parameter, private :: f8 ='(''tbt: '',a,t53,''='',f10.4)'
  character(len=*), parameter, private :: f9 ='(''tbt: '',a,t53,''='',tr1,e9.3)'
  character(len=*), parameter, private :: f15='(''tbt: '',a,t53,''='',2(tr1,i0,'' x''),'' '',i0)'


contains

  subroutine read_tbt_generic(na_u, lasto)

    use fdf, only: fdf_defined
    use m_ts_method, only: ts_init_regions
    use m_tbt_diag, only: init_diag

    ! The number of atoms
    integer, intent(in) :: na_u
    ! A summated list of last orbitals on atoms.
    integer, intent(in) :: lasto(0:na_u)

    ! Initialize the buffer regions
    if ( fdf_defined('TBT.Atoms.Buffer') ) then
       call ts_init_regions('TBT',na_u,lasto)
    else
       call ts_init_regions('TS',na_u,lasto)
    end if

    ! Initialize the diagonalization method.
    call init_diag( )

  end subroutine read_tbt_generic


  ! > Reads the chemical potentials as well as the applied
  ! Bias.
  ! The bias is an intricate part of the chemical potential why it
  ! is read in here.
  subroutine read_tbt_chem_pot( )

    use fdf, only : fdf_get
    use units, only: eV, Kelvin

    use m_ts_chem_pot, only : fdf_nmu, fdffake_mu, fdf_mu, name
    
    use m_tbt_hs, only: Volt

    implicit none

    ! *******************
    ! * LOCAL variables *
    ! *******************
    logical :: err
    integer :: i

    ! Read in the temperature
    kT = fdf_get('ElectronicTemperature',1.9e-3_dp,'Ry')
    kT = fdf_get('TS.ElectronicTemperature',kT,'Ry')
    kT = fdf_get('TBT.ElectronicTemperature',kT,'Ry')

    ! Read in the chemical potentials
    N_mu = fdf_nmu('TBT',kT,mus)
    if ( N_mu < 1 ) then
       N_mu = fdf_nmu('TS',kT,mus)
    end if
    err = .true.
    if ( N_mu < 1 ) then
       err = .false.
       N_mu = fdffake_mu(mus,kT,Volt)
    end if
    do i = 1 , N_mu
       ! Default things that could be of importance
       if ( fdf_mu('TBT',mus(i),kT,Volt) ) then
          ! success
       else if ( fdf_mu('TS',mus(i),kT,Volt) ) then
          ! success
       else if ( err ) then
          ! only error out if it couldn't be found and forced
          ! created
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
               &potentials. I.e. you cannot lift the frequency spectra.')
       end if
    end do
#endif

  end subroutine read_tbt_chem_pot


  ! Reads all information regarding the electrodes, nothing more.
  subroutine read_tbt_elec( cell, na_u, xa, lasto)

    use fdf, only : fdf_get, fdf_obsolete, fdf_deprecated, leqi
    use parallel, only : IONode
    use intrinsic_missing, only : IDX_SPC_PROJ, EYE

    use m_os, only : file_exist

    use files, only: slabel
    use units, only: eV

    use m_tbt_hs, only: spin_idx

    use m_ts_chem_pot, only : copy, chem_pot_add_Elec

    use m_ts_electype, only : fdf_nElec, fdf_Elec
    use m_ts_electype, only : Name, TotUsedOrbs, TotUsedAtoms
    use m_ts_electype, only : init_Elec_sim

    use m_ts_method, only : ts_init_electrodes, a_isBuffer

    implicit none

    ! *******************
    ! * INPUT variables *
    ! *******************
    real(dp), intent(in) :: cell(3,3)
    integer,  intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)

    ! *******************
    ! * LOCAL variables *
    ! *******************
    integer :: i, j
    real(dp) :: tmp33(3,3), rtmp
    logical :: err
    character(len=200) :: chars

    if ( N_mu == 0 ) call die('read_tbt_elecs: error in programming')

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
    Elecs(:)%Bulk = fdf_get('TS.Elecs.Bulk',.true.) ! default everything to bulk electrodes
    Elecs(:)%Bulk = fdf_get('TBT.Elecs.Bulk',Elecs(1)%Bulk)

    rtmp = fdf_get('TS.Elecs.Eta',0.0001*eV,'Ry')
    rtmp = fdf_get('TBT.Elecs.Eta',rtmp,'Ry')
    Elecs(:)%Eta = rtmp

    ! whether all calculations should be performed
    ! "out-of-core" i.e. whether the GF files should be created or not
    ! In tbtrans this is now defaulted to in-core
    Elecs(:)%out_of_core = fdf_get('TBT.Elecs.Out-of-core',.false.)

    ! Whether we should try and re-use the surface Green function 
    ! files
    Elecs(:)%ReUseGF = fdf_get('TS.Elecs.GF.ReUse',.true.)
    Elecs(:)%ReUseGF = fdf_get('TBT.Elecs.GF.ReUse',Elecs(1)%ReUseGF)

    ! Will stop after creating the GF files.
    stop_after_GS = fdf_get('TBT.Elecs.GF.Only',.false.)

    do i = 1 , N_Elec

       ! If we only have 2 electrodes we take them 
       ! as though the atomic indices are the first and last
       ! respectively.
       if ( N_Elec == 2 ) then
          if ( i == 1 ) then
             err = fdf_Elec('TBT',slabel,Elecs(i),N_mu,mus,idx_a= 1, &
                  name_prefix = name_prefix)
             if ( .not. err ) &
                  err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus,idx_a= 1, &
                  name_prefix = name_prefix)
          else
             err = fdf_Elec('TBT',slabel,Elecs(i),N_mu,mus,idx_a=-1, &
                  name_prefix = name_prefix)
             if ( .not. err ) &
                  err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus,idx_a=-1, &
                  name_prefix = name_prefix)
          end if
       else
          ! Default things that could be of importance
          err = fdf_Elec('TBT',slabel,Elecs(i),N_mu,mus)
          if ( .not. err ) &
               err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus, &
               name_prefix = name_prefix)
       end if
       if ( .not. err ) then
          call die('Could not find electrode: '//trim(name(Elecs(i))))
       end if
       if ( Elecs(i)%idx_a < 0 ) &
            Elecs(i)%idx_a = na_u + Elecs(i)%idx_a + 1
       if ( Elecs(i)%idx_a < 1 .or. &
            na_u < Elecs(i)%idx_a ) then
          print *,Elecs(i)%idx_a,na_u
          call die("Electrode position does not exist")
       end if
       if ( N_Elec == 2 ) then
          ! Correct for buffer atoms, first electrode steps "up"
          ! second electrode steps "down"
          if ( i == 1 ) then
             j = Elecs(i)%idx_a
             do while ( a_isBuffer(j) )
                j = j + 1
             end do
             Elecs(i)%idx_a = j
          else
             j = Elecs(i)%idx_a + TotUsedAtoms(Elecs(i)) - 1
             do while ( a_isBuffer(j) )
                j = j - 1
             end do
             Elecs(i)%idx_a = j - TotUsedAtoms(Elecs(i)) + 1
          end if
       end if
       ! set the placement in orbitals
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

       ! Initialize electrode parameters
       call init_Elec_sim(Elecs(i),cell,na_u,xa)

    end do

    ! Initialize the electrode regions
    call ts_init_electrodes(na_u,lasto,N_Elec,Elecs)

    ! If many electrodes, no transport direction can be specified
    ! Hence we use this as an error-check (also for N_Elec == 1)
    if ( N_Elec /= 2 ) then
       ! Signals no specific unit-cell direction of transport
       ts_tdir = - N_Elec
       ts_tidx = - N_Elec
    else

       ! Retrieve the indices of the unit-cell directions
       ! according to the electrode transport directions.
       ! We have already calculated the pivoting table for
       ! the electrodes
       i = Elecs(1)%pvt(Elecs(1)%t_dir)
       j = Elecs(2)%pvt(Elecs(2)%t_dir)

       if ( i == j ) then
          ! The transport direction for the electrodes are the same...
          ts_tidx = i

          ! Calculate Cartesian transport direction
          call eye(3,tmp33)
          ts_tdir = IDX_SPC_PROJ(tmp33,cell(:,ts_tidx))

       else
          ! In case we have a skewed transport direction
          ! we have some restrictions...
          ts_tidx = - N_Elec
          ts_tdir = - N_Elec
       end if

    end if

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
       write(*,'(a)') 'Please stop this madness. What where you thinking?'
       call die('Electrodes occupy the entire device!!!')
    end if

  end subroutine read_tbt_elec


  subroutine read_tbt_after_Elec(nspin, cell, na_u, lasto, xa, kscell, kdispl)

    use fdf, only: fdf_get, leqi
    use parallel, only: IONode

    use m_ts_method, only: ts_method, TS_BTD

    use m_tbt_contour, only: read_contour_options

    use m_tbt_save, only: init_save_options
    use m_tbt_sigma_save, only: init_Sigma_options
    use m_tbt_dH, only: init_dh_options
    
    ! *******************
    ! * INPUT variables *
    ! *******************
    integer, intent(in) :: nspin
    real(dp), intent(in) :: cell(3,3)
    integer,  intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer,  intent(in) :: kscell(3,3)
    real(dp), intent(in) :: kdispl(3)

    integer :: i
    logical :: ltmp, Gamma3(3)
    character(len=150) :: chars
    
    ! we must have read the electrodes first
    if ( N_Elec == 0 ) call die('read_tbt_options: Error in programming')

    percent_tracker = fdf_get('TBT.Progress',5)
    percent_tracker = max(1,percent_tracker)

    ! Reading the Transiesta solution method
    chars = fdf_get('TBT.SolutionMethod','BTD')
    if ( leqi(chars,'BTD').or.leqi(chars,'tri') ) then
       ts_method = TS_BTD
    else
       call die('Unrecognized TBtrans solution method: '//trim(chars))
    end if

    chars = fdf_get('TS.BTD.Optimize','speed')
    chars = fdf_get('TBT.BTD.Optimize',trim(chars))
    if ( leqi(chars,'speed') ) then
       BTD_method = 0
    else if ( leqi(chars,'memory') ) then
       BTD_method = 1
    else
       call die('Could not determine flag TBT.BTD.Optimize, please &
            &see manual.')
    end if

    ! Whether we should assert and calculate
    ! all transmission amplitudes
    ltmp = fdf_get('TBT.T.Elecs.All',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('T-all'.kv.1)
    end if

    N_eigen = fdf_get('TBT.T.Eig',0)
    if ( N_eigen > 0 ) then
       save_DATA = save_DATA // ('T-eig'.kv.N_eigen)
    end if

    ! Should we calculate DOS of electrode bulk Green function
    ltmp = fdf_get('TBT.DOS.Elecs', .false. )
    if ( ltmp ) then
       save_DATA = save_DATA // ('DOS-Elecs'.kv.1)
    end if

    ! Should we calculate DOS of Green function
    ltmp = fdf_get('TBT.DOS.Gf', N_Elec == 1 )
    if ( ltmp ) then
       save_DATA = save_DATA // ('DOS-Gf'.kv.1)
    end if

    ! Should we calculate DOS of spectral function
    ltmp = fdf_get('TBT.DOS.A',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('DOS-A'.kv.1)
    end if

    ! Should we calculate DOS of all spectral functions
    ltmp = fdf_get('TBT.DOS.A.All',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('DOS-A-all'.kv.1)
       save_DATA = save_DATA // ('DOS-A'.kv.1)
    end if

    ! Should we calculate orbital current
    ltmp = fdf_get('TBT.Current.Orb', .false. )
    if ( ltmp .and. ('DOS-A'.in.save_DATA)) then
       save_DATA = save_DATA // ('orb-current'.kv.1)
    else if ( ltmp .and. IONode ) then
       write(*,'(a,/,a)')'WARNING: Will not calculate the orbital currents, &
            &the spectral function needs to be calculated for this to &
            &apply.','Set TBT.DOS.A T to calculate orbital currents.'
    end if

    ltmp = fdf_get('TBT.T.Out',N_Elec == 1)
    if ( ltmp ) then
       save_DATA = save_DATA // ('T-sum-out'.kv.1)
    end if

#ifdef NCDF_4
    call init_Sigma_options( save_DATA )

    call init_dH_options( )
#endif
    call init_save_options( )

    ! read in contour options
    call read_contour_options( N_Elec, Elecs, N_mu, mus )

    ! Check for Gamma in each direction
    do i = 1 , 3
       if ( kdispl(i) /= 0._dp ) then
          Gamma3(i) = .false.
       else if ( sum(kscell(i,:)) > 1 ) then
          ! Note it is the off-diagonal for this unit-cell
          ! direction
          Gamma3(i) = .false.
       else
          Gamma3(i) = .true.
       end if
    end do

    do i = 1 , N_Elec
       ! Initialize the electrode quantities for the stored values
       call check_Elec_sim(Elecs(i), nspin, cell, na_u, xa, &
            Elecs_xa_EPS, lasto, Gamma3)
    end do

  end subroutine read_tbt_after_Elec

  subroutine print_tbt_options(nspin)

    use units, only: Kelvin, eV
    use parallel, only: IONode
    use files, only: slabel

    use m_tbt_contour, only: print_contour_tbt_options, io_contour_tbt
    use m_tbt_contour, only: print_contour_tbt_block
    use m_tbt_save, only: print_save_options
    use m_tbt_diag, only: print_diag
#ifdef NCDF_4
    use m_tbt_dH, only: print_dH_options
    use m_tbt_sigma_save, only: print_Sigma_options
#endif
    use m_tbt_hs, only: Volt, IsVolt, spin_idx

    integer, intent(in) :: nspin
    
    integer :: i
    
    if ( .not. IONode ) return

    write(*,*)
    write(*,f11) repeat('*', 62)

    write(*,f6) 'Electronic temperature',kT/Kelvin,'K'
    if ( IsVolt ) then
       write(*,f6) 'Voltage', Volt/eV,'Volts'
    else
       write(*,f11) 'No applied bias'
    end if
    write(*,f1) 'Saving DOS from bulk electrodes',('DOS-Elecs'.in.save_DATA)
    write(*,f1) 'Saving DOS from Green function',('DOS-Gf'.in.save_DATA)
    if ( 'DOS-A-all' .in. save_DATA ) then
       write(*,f1) 'Saving DOS from all spectral functions',.true.
    else
       write(*,f1) 'Saving DOS from spectral functions',('DOS-A' .in. save_DATA)
    end if
    write(*,f12) 'Calc. # transmission eigenvalues',N_eigen
    write(*,f1) 'Calc. T between all electrodes',('T-all'.in.save_DATA)
    write(*,f1) 'Calc. total T out of electrodes',('T-sum-out'.in.save_DATA)
    if ( nspin > 1 ) then
       if ( spin_idx == 0 ) then
          write(*,f11) 'Calculate all spin-channels'
       else
          write(*,f9) 'Calculate spin-channel',spin_idx
       end if
    else
       write(*,f11) 'Single spin Hamiltonian'
    end if

    call print_diag()
#ifdef NCDF_4
    call print_Sigma_options( save_DATA )
    call print_dH_options()
#endif
    call print_save_options()

    write(*,f11)'          >> Electrodes << '
    do i = 1 , size(Elecs)
       call print_settings(Elecs(i),'tbt')
    end do
    
    call print_contour_tbt_options( 'TBT' )
    
    write(*,f11) repeat('*', 62)
    write(*,*)
    
    call io_contour_tbt(slabel)

    write(*,f11) repeat('<', 62)

    call print_contour_tbt_block( 'TBT' )

    write(*,f11) repeat('<', 62)

  end subroutine print_tbt_options

  subroutine print_tbt_warnings( Gamma )

    use fdf, only: fdf_get
    use units, only: eV
    use parallel, only: IONode

    use m_tbt_hs, only: Volt

    ! Whether the user requests a Gamma calculation
    logical, intent(in) :: Gamma

    integer :: i
    logical :: ltmp

    if ( N_eigen < 0 ) then
       call die('Number of transmission eigenvalues MUST be &
            &zero or positive.')
    end if

    if ( ('orb-current' .in.save_DATA) ) then
       ltmp = .not. fdf_get('SpinSpiral',.false.)
       ltmp = fdf_get('TBT.Symmetry.TimeReversal',ltmp)
       if ( IONode .and. .not. Gamma ) then
          write(*,'(a,/,a)') 'WARNING: k-averaging orbital currents with &
               &time-reversal symmetry will not reproduce','the correct &
               &orbital current. Set TBT.Symmetry.TimeReversal F'
       end if
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

  end subroutine print_tbt_warnings

end module m_tbt_options
