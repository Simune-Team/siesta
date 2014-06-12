module m_ts_options

  use precision, only : dp
  use siesta_options, only : FixSpin, isolve, SOLVE_TRANSI
  use sys, only : die

  use m_ts_electype
  use m_ts_chem_pot
  use m_ts_tdir

  implicit none

  public
  save

  ! Controls to save the TSHS file
  logical  :: SaveTSHS = .true. 
  ! whether we should only save the overlap matricx
  logical  :: onlyS = .false. 
  ! whether we will use the bias-contour
  logical  :: IsVolt = .false.
  ! maximum difference between chemical potentials
  real(dp) :: Volt = 0._dp
  ! Electrodes and different chemical potentials
  integer :: N_Elec = 0
  type(Elec), allocatable, target :: Elecs(:)
  integer :: N_mu = 0
  type(ts_mu), allocatable, target :: mus(:)
  
  logical :: ImmediateTSmode = .false. ! will determine to immediately start the transiesta
                                       ! SCF. This is useful when you already have a converged
                                       ! siesta DM

  ! Whether we should remove the inner-cell distances
  logical :: RemUCellDistance = .false.
  ! Flag to control whether we should update the forces (i.e. calculate energy-density matrix)
  logical :: Calc_Forces = .true.

  ! If the energy-contour is not perfectly divisable by the number of nodes then adjust
  integer :: opt_TriMat_method = 0 ! Optimization method for determining the best tri-diagonal matrix split
  ! 0  == We optimize for speed
  ! 1  == We optimize for memory

  ! Determines whether the voltage-drop should be located in the constriction
  ! I.e. if the electrode starts at 10 Ang and the central region ends at 20 Ang
  ! then the voltage drop will only take place between 10.125 Ang and 19.875 Ang
  logical :: VoltageInC = .false.

  ! A quantity describing the accuracy of the coordinates of the 
  ! electrodes.
  ! * Should only be edited by experienced users *
  real(dp) :: Elecs_xa_EPS = 1.e-4_dp

  ! The mixing weight in the transiesta cycles...
  real(dp) :: ts_wmix ! = wmix

  ! The user can request to analyze the system, returning information about the 
  ! tri-diagonalization partition and the contour
  logical :: TS_Analyze = .false.
  integer :: TS_bandwidth_algo = 0
  
contains
  
  subroutine read_ts_options( wmix, kT, ucell, na_u, xa, lasto)

    use alloc
    use files, only : slabel
    use fdf, only : fdf_get, fdf_deprecated, fdf_obsolete
    use fdf, only : leqi
    use parallel, only: IOnode, Nodes
    use units, only: eV, Ang, Kelvin

    use m_ts_cctype
    use m_ts_global_vars, only : TSmode, ts_istep
    use m_ts_io, only : ts_read_TSHS_opt

    use m_ts_contour
    use m_ts_contour_eq,  only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E

    use m_ts_io_contour
    use m_ts_method
    use m_ts_weight
    use m_ts_charge
    use m_ts_tdir

#ifdef MUMPS
    use m_ts_mumps_init, only : MUMPS_mem, MUMPS_ordering, MUMPS_block
#endif
    
    use m_monitor
    use m_bandwidth

    implicit none
    
! *******************
! * INPUT variables *
! *******************
    real(dp), intent(in) :: wmix, kT
    real(dp),intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)

! *******************
! * LOCAL variables *
! *******************
    real(dp) :: tmp
    logical :: err
    character(len=200) :: c, chars
    integer :: i, j, idx, idx1, idx2

    type(ts_mu) :: tmp_mu
    ! External routines
    real(dp) :: dot
    external :: dot

    if (IOnode) then
       write(*,*)
       write(*,11) repeat('*', 62)
    end if

    ! Set ts_istep default
    ts_istep = 0

    ! Read in general values that should be used in the electrode generation
    ! I.e. these FDF-parameters are used for diagon runs with transiesta
    saveTSHS = fdf_get('TS.SaveHS',.true.)
    onlyS    = fdf_get('TS.onlyS',.false.)

    if ( SaveTSHS .and. FixSpin ) then
       write(*,*) 'Fixed spin not possible with Transiesta!'
       write(*,*) 'Electrodes with fixed spin is not possible with Transiesta !'
       call die('Stopping code')
    end if

    if ( .not. TSmode ) then
       if ( IONode ) then
          write(*,1) 'Save H and S matrices', saveTSHS
          write(*,1) 'Save S and quit (onlyS)', onlyS
          write(*,11) repeat('*', 62)
          write(*,*)
       end if
       return
    end if

    ! Read in the mixing for the transiesta cycles
    ts_wmix = fdf_get('TS.MixingWeight',wmix)
    
    ! Read in the transport direction
    chars = fdf_get('TS.TransportDirection','c')
    if ( leqi(chars,'a') .or. leqi(chars,'a1') ) then
       ts_tdir = 1
    else if ( leqi(chars,'b') .or. leqi(chars,'a2') ) then
       ts_tdir = 2
    else if ( leqi(chars,'c') .or. leqi(chars,'a3') ) then
       ts_tdir = 3
    else
       call die('Transport direction not in [a|b|c|a1|a2|a3]')
    end if

    ! Read in information about the voltage placement.
    chars = fdf_get('TS.HartreePotential.Position','central')
    VoltageInC = .true.
    if ( leqi(trim(chars),'cell') ) then
       VoltageInC = .false.
    else if ( leqi(trim(chars),'central') .or. &
         leqi(trim(chars),'scat') ) then
       VoltageInC = .true.
    end if

    ! Reading the Transiesta solution method
    chars = fdf_get('TS.SolutionMethod','tri')
    if ( leqi(chars,'sparse') ) then
       ts_method = TS_SPARSITY
    else if ( leqi(chars,'tri') ) then
       ts_method = TS_SPARSITY_TRI
#ifdef MUMPS
    else if ( leqi(chars,'mumps') ) then
       ts_method = TS_SPARSITY_MUMPS
#endif
    else
       call die('Unrecognized Transiesta solution method: '//trim(chars))
    end if

#ifdef MUMPS
    MUMPS_mem   = fdf_get('TS.MUMPS.Mem',20)
    MUMPS_block = fdf_get('TS.MUMPS.BlockingFactor',-8)
    chars = fdf_get('TS.MUMPS.Ordering','auto')
    if ( leqi(chars,'auto') ) then
       MUMPS_ordering = 7
    else if ( leqi(chars,'amd') ) then
       MUMPS_ordering = 0
    else if ( leqi(chars,'amf') ) then
       MUMPS_ordering = 2
    else if ( leqi(chars,'scotch') ) then
       MUMPS_ordering = 3
    else if ( leqi(chars,'pord') ) then
       MUMPS_ordering = 4
    else if ( leqi(chars,'metis') ) then
       MUMPS_ordering = 5
    else if ( leqi(chars,'qamd') ) then
       MUMPS_ordering = 6
    else
       call die('Unknown MUMPS ordering.')
    end if
#endif


    ! currently this does not work
    !ImmediateTSmode = fdf_get('TS.SCFImmediate',.false.)

    chars = fdf_get('TS.TriMat.Optimize','speed')
    if ( leqi(chars,'speed') ) then
       opt_TriMat_method = 0
    else if ( leqi(chars,'memory') ) then
       opt_TriMat_method = 1
    else
       call die('Could not determine flag TS.TriMat.Optimize, please &
            &see manual.')
    end if

    ! Determine whether the user wishes to only do an analyzation
    TS_Analyze = fdf_get('TS.Analyze',.false.)
    if ( TS_Analyze ) then
       ! Default
       TS_bandwidth_algo = 0
       chars = fdf_get('TS.BandwidthReduction','reverse-CutHill-Mckee')
       i = index(chars,'-')
       if ( i > 0 ) then
          if ( leqi(chars(1:i-1),'reverse') .or. &
               leqi(chars(1:i-1),'rev') ) then
             TS_bandwidth_algo = BW_REVERSE
             chars = chars(i+1:)
          end if
       end if
       if ( leqi(chars,'cuthill-mckee') .or. &
            leqi(chars,'cm') ) then
          TS_bandwidth_algo = TS_bandwidth_algo + BW_CUTHILL_MCKEE
       else if ( leqi(chars,'cuthill-mckee-priority') .or. &
            leqi(chars,'cm-p') ) then
          TS_bandwidth_algo = TS_bandwidth_algo + BW_CUTHILL_MCKEE_PRIORITY
       else if ( leqi(chars,'papior') ) then
          TS_bandwidth_algo = TS_bandwidth_algo + BW_PAPIOR
       else
          call die('Unrecognized option for Bandwidth algorithm: '//trim(chars))
       end if
    end if
    
    chars = fdf_get('TS.ChargeCorrection','none')
    TS_RHOCORR_METHOD = 0
    if ( leqi(chars,'none') ) then
       TS_RHOCORR_METHOD = 0
    else if ( leqi(chars,'b') .or. leqi(chars,'buffer') ) then
       TS_RHOCORR_METHOD = TS_RHOCORR_BUFFER
    !else if ( leqi(chars,'up') .or. leqi(chars,'update') ) then
    !   TS_RHOCORR_METHOD = TS_RHOCORR_UPDATE
    end if
    TS_RHOCORR_FACTOR = fdf_get('TS.ChargeCorrection.Factor',.75_dp)
    if ( TS_RHOCORR_FACTOR < 0.0_dp .or. &
         1.0_dp < TS_RHOCORR_FACTOR) then
       call die("Charge correction factor must be in the range [0;1]")
    endif
    
    ! whether to calculate the forces or not (default calculate everything)
    Calc_Forces = fdf_get('TS.Forces',.true.)

    ! The sign can not be chosen from this (several mu, where to define it)
    Volt   = fdf_get('TS.Voltage',0._dp,'Ry') 
    ! Voltage situation is above 0.01 mV
    IsVolt = abs(Volt/eV) > 0.00001_dp

    ! This should never be used!!!!
    ! It is required to fix the potential in the cell
    call fdf_obsolete('TS.UseVFix')

    ! { Entering electrode land... First we describe all the deprecated options

    ! the title of the green's functions are now non-generic
    call fdf_obsolete('TS.GFTitle')


    ! The regular options for describing the electrodes can not be 
    ! used anymore...
    do i = 1 , 2
       if ( i == 1 ) then
          chars = 'Left'
       else
          chars = 'Right'
       end if
       call fdf_obsolete('TS.HSFile'//trim(chars))
       call fdf_obsolete('TS.GFFile'//trim(chars))
       call fdf_obsolete('TS.NumUsedAtoms'//trim(chars))
       call fdf_obsolete('TS.ReplicateA1'//trim(chars))
       call fdf_obsolete('TS.ReplicateA2'//trim(chars))
    end do

    ! notice that this does not have the same meaning... 
    call fdf_deprecated('TS.UpdateDMCROnly','TS.Elecs.DM.CrossTerms')
    call fdf_deprecated('TS.UseBulk','TS.Elecs.Bulk')

    ! Read in the chemical potentials
    N_mu = fdf_nmu('TS',mus)
    if ( N_mu < 1 ) then
       N_mu = fdffake_mu(mus,Volt)
    else
       do i = 1 , N_mu
          ! Default things that could be of importance
          if ( .not. fdf_mu('TS',mus(i),Volt) ) then
             call die('Could not find chemical potential: ' &
                  //trim(name(mus(i))))
          end if
          ! Attach the ID
          mus(i)%ID = i
       end do
    end if

    ! To determine the same coordinate nature of the electrodes
    Elecs_xa_EPS= fdf_get('TS.Elecs.Coord.Eps',1.e-4_dp,'Bohr')

    ! detect how many electrodes we have
    N_Elec = fdf_nElec('TS',Elecs)
    if ( N_Elec < 1 ) then
       call die('Please see the manual for how to construct an example &
            &electrode configuration (or use Util/TS/tselecs.sh)')
    end if
    ! If only one electrode you are not allowed to move the fermi-level
    ! of the electrode. That should be done by other means (i.e. use NetCharge)
    if ( N_Elec == 1 ) then
       ! Notice that below the chemical potential gets corrected
       ! EVEN if the user supplied a bias.
       if ( IsVolt .and. IONode ) then
          c = '(''transiesta: ***'',a)'
          write(*,c) 'Single electrode calculations does not allow shifting the chemical potential.'
          write(*,c) 'You should do that by changing the states filled in the system.'
          write(*,c) 'Consult the manual of how to do this.'
       end if
       IsVolt = .false.
    end if

    ! If many electrodes, no transport direction can be specified
    ! Hence we use this as an error-check (also for N_Elec == 1)
    if ( N_Elec /= 2 ) then
       ts_tdir = - N_Elec
    end if

    ! Setup default parameters for the electrodes
    ! first electrode is the "left"
    ! last electrode is the "right"
    ! the remaining electrodes have their chemical potential at 0
    ! Currently the transport direction for all electrodes is the default
    ! We should probably warn if +2 electrodes are used and t_dir is the
    ! same for all electrodes... Then the user needs to know what (s)he is doing...
    Elecs(:)%t_dir = ts_tdir
    Elecs(:)%Bulk  = fdf_get('TS.Elecs.Bulk',.true.) ! default everything to bulk electrodes
    if ( .not. Elecs(1)%Bulk ) then
       Elecs(:)%DM_CrossTerms = .true.
    else
       ! default to not update the cross-terms
       Elecs(:)%DM_CrossTerms = fdf_get('TS.Elecs.DM.CrossTerms',.false.)
    end if
    ! We default to not calculate the band-bottom...
    ! TODO move to TS.Analyze step..., no need to have this in TS-scheme...
    Elecs(:)%BandBottom = fdf_get('TS.Elecs.BandBottom', .false.)
    ! whether or not the electrodes should be re-instantiated
    call fdf_deprecated('TS.CalcGF','TS.Elecs.GF.ReUse')
    call fdf_deprecated('TS.ReUseGF','TS.Elecs.GF.ReUse')
    err = fdf_get('TS.ReUseGF',.false.)
    Elecs(:)%ReUseGF = fdf_get('TS.Elecs.GF.ReUse',err)

    ! whether all calculations should be performed
    ! "out-of-core" i.e. whether the GF files should be created or not
    Elecs(:)%out_of_core = fdf_get('TS.Elecs.Out-of-core',.true.)

    do i = 1 , N_Elec
       ! Default things that could be of importance
       if ( .not. fdf_Elec('TS',slabel,Elecs(i),N_mu,mus) ) then
          call die('Could not find electrode: '//trim(name(Elecs(i))))
       end if
       ! set the placement in orbitals
       if ( Elecs(i)%idx_na < 0 ) &
            Elecs(i)%idx_na = na_u + Elecs(i)%idx_na + 1
       if ( Elecs(i)%idx_na < 1 .or. &
            na_u < Elecs(i)%idx_na ) &
            call die("Electrode position does not exist")
       Elecs(i)%idx_no = lasto(Elecs(i)%idx_na-1)+1

    end do

    ! Check that the current transport direction is "aligned"
    ! TODO when we can deal with arbitrary electrodes this should be altered
    if ( N_Elec == 2 .and. &
         all(Elecs(1)%t_dir == Elecs(:)%t_dir) ) then
       ts_tdir = Elecs(1)%t_dir
    end if

    if ( .not. IsVolt ) then
       ! force it to be zero... can be necessary if considering single electrode
       ! calculations (assures V == 0)
       Volt = 0._dp

       ! copy over electrode...
       call copy(mus(1),tmp_mu)
       
       ! Deallocate all strings
       do i = 1 , N_mu
          deallocate(mus(i)%Eq_seg)
       end do
       deallocate(mus)

       ! create the first chemical potential again
       N_mu = 1
       allocate(mus(1))
       call copy(tmp_mu,mus(1))
       deallocate(tmp_mu%Eq_seg)

       ! Firmly assure the chemical potential to be zero
       mus(1)%mu = 0._dp
       mus(1)%ID = 1
       
       ! Assign all electrodes to the same chemical potential
       do i = 1 , N_Elec
          Elecs(i)%mu => mus(1)
       end do

    end if

    ! We can now check whether two chemical potentials are the same
    ! We must do this after checking for the equilibrium case as we wish to
    ! allow users to retain all chemical potentials at 0 eV
    do i = 1 , N_mu - 1
       do j = i + 1 , N_mu
          if ( abs(mus(i)%mu - mus(j)%mu) < 0.00001_dp*eV ) then
             call die('Two chemical potentials: '//trim(name(mus(i)))//' and ' &
                  //trim(name(mus(j)))//' are the same, in bias calculations this &
                  &is not allowed.')
          end if
       end do
    end do

    ! Populate the electrodes in the chemical potential type
    do i = 1 , N_Elec
       do j = 1 , N_mu
          if ( associated(Elecs(i)%mu,target=mus(j)) ) then
             call chem_pot_add_Elec(mus(j),i)
             exit
          end if
       end do
    end do

    ! check that all electrodes and chemical potentials are paired in
    ! some way.
    if ( any(mus(:)%N_El == 0) ) then
       call die('A/Some chemical potential(s) has not been assigned any electrodes. &
            &All chemical potentials *MUST* be assigned an electrode')
    end if

    ! check that all have at least 2 contour points on the equilibrium contour
    ! the 3rd is the fictive pole segment
    if ( .not. all(Eq_segs(mus(:)) > 2) ) then
       call die('All chemical potentials does not have at least 2 equilibrium contours')
    end if
    
    if ( na_u <= sum(TotUsedAtoms(Elecs)) ) then
       call die('Electrodes occupy the entire device!!!')
    end if

    err = .false.
    ! we need to check that they indeed do not overlap
    do i = 1 , N_Elec
       idx1 = Elecs(i)%idx_na
       idx2 = idx1 + TotUsedAtoms(Elecs(i)) - 1
       ! we need to check every electrode,
       ! specifically because if one of the electrodes is fully located
       ! inside the other and we check the "small" one 
       do j = 1 , N_Elec
          if ( i == j ) cycle
          idx = Elecs(j)%idx_na
          if ( (idx <= idx1 .and. &
               idx1 < idx + TotUsedAtoms(Elecs(j))) ) then
             err = .true.
          end if
          if ( (idx <= idx2 .and. &
               idx2 < idx + TotUsedAtoms(Elecs(j))) ) then
             err = .true.
          end if
          if ( err ) then
             write(*,*) 'Electrode: '//trim(Name(Elecs(i)))
             write(*,'(a,i0,a,i0)') 'Positions: ',idx1,' -- ',idx2 
             idx1 = Elecs(j)%idx_na
             idx2 = idx1 + TotUsedAtoms(Elecs(j)) - 1
             write(*,*) 'Electrode: '//trim(Name(Elecs(j)))
             write(*,'(a,i0,a,i0)') 'Positions: ',idx1,' -- ',idx2 
             call die('Overlapping electrodes is not physical, please correct.')
          end if
       end do
    end do

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

    ! WILL WORK EVENTUALLY
    if ( N_Elec > 2 .and. IsVolt ) call die('Several electrodes and bias does not work')

    ! Update the weight function
    chars = fdf_get('TS.Weight.k.Method','correlated')
    if ( leqi(chars,'correlated') ) then
       TS_W_K_METHOD = TS_W_K_CORRELATED
    else if ( leqi(chars,'half-correlated') ) then
       TS_W_K_METHOD = TS_W_K_HALF_CORRELATED
       call die('Currently not functioning')
    else if ( leqi(chars,'uncorrelated') ) then
       TS_W_K_METHOD = TS_W_K_UNCORRELATED
    else
       call die('Could not determine flag TS.Weight.k.Method, &
            &please see manual.')
    end if

    ! The default weighting method is correlated if
    ! atom-atom is utilised
    TS_W_METHOD = TS_W_CORRELATED
    chars = fdf_get('TS.Weight.Method','orb-orb')
    ! first check whether we have correlated weighting
    i = index(chars,'+')
    if ( i > 0 ) then
       ! we do have something else
       if ( leqi(chars(1:i-1),'correlated') .or. &
            leqi(chars(1:i-1),'corr') ) then
          TS_W_METHOD = TS_W_CORRELATED
       else if ( leqi(chars(1:i-1),'uncorrelated') .or. &
            leqi(chars(1:i-1),'uncorr') ) then
          TS_W_METHOD = 0 ! non-correlated
       else
          call die('Unrecognized second option for TS.Weight.Method &
               &must be [[un]correlated+][orb-orb|tr-atom-atom|sum-atom-atom]')
       end if
       chars = chars(i+1:)
    end if
    if ( leqi(chars,'orb-orb') ) then
       TS_W_METHOD = TS_W_ORB_ORB
       ! this does not make sense to make correlated, hence always assign
    else if ( leqi(chars,'tr-atom-atom') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_TR_ATOM_ATOM 
    else if ( leqi(chars,'tr-atom-orb') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_TR_ATOM_ORB
    else if ( leqi(chars,'sum-atom-atom') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_SUM_ATOM_ATOM 
    else if ( leqi(chars,'sum-atom-orb') ) then
       TS_W_METHOD = TS_W_METHOD + TS_W_SUM_ATOM_ORB
    else
       call die('Unrecognized option for TS.Weight.Method &
            &must be [[un]correlated+|][orb-orb|tr-atom-[atom|orb]|sum-atom-[atom|orb]]')
    end if
    if ( TS_W_METHOD /= TS_W_ORB_ORB ) then
       ! We do not allow to do the half-correlated,
       if ( TS_W_K_METHOD == TS_W_K_HALF_CORRELATED ) then
          call die('The uncorrelated weighting does not work &
               &with trace-weighting of the density matrix.')
       end if
       
    end if


    ! read in contour options
    if ( TSmode ) then
       call read_contour_options( N_Elec, Elecs, N_mu, mus, kT, IsVolt, Volt )
    end if

    ! read in buffer information
    if ( TSMode ) then
       ! initialize regions of the electrodes and device
       ! the number of LCAO orbitals on each atom will not change
       call ts_init_regions('TS',N_Elec,Elecs,na_u,lasto)
    end if

    ! Show the deprecated and obsolete labels
    call fdf_deprecated('TS.TriDiag','TS.SolutionMethod')
    call fdf_obsolete('TS.FixContactCharge')
    call fdf_obsolete('TS.KxyPoints')
    call fdf_obsolete('TS.NKVoltScale')

    if (IONode .and. TSmode ) then
       write(*,1) 'Save H and S matrices', saveTSHS
       if ( TS_Analyze ) then
          write(*,11)'Will analyze bandwidth of LCAO sparse matrix and quit'
          chars = ''
          i = TS_bandwidth_algo
          if ( TS_bandwidth_algo >= BW_REVERSE ) then
             chars = 'reverse-'
             i = TS_bandwidth_algo - BW_REVERSE
          end if
          select case ( i ) 
          case ( BW_CUTHILL_MCKEE ) 
             chars = trim(chars)//'Cuthill-Mckee'
          case ( BW_CUTHILL_MCKEE_PRIORITY ) 
             chars = trim(chars)//'Cuthill-Mckee with T-dir priority'
          case ( BW_PAPIOR ) 
             chars = trim(chars)//'Papior'
          case default 
             call die('Unrecognized bandwidth algorithm')
          end select
          write(*,10)'Bandwidth algorithm',trim(chars)
       end if
       write(*,7) 'Electronic temperature',kT/Kelvin,'K'
       write(chars,'(a,i0)') 'A',ts_tdir
       write(*,10) 'Transport along unit-cell vector',trim(chars)
       if ( ts_method == TS_SPARSITY ) then
          write(*,10)'Solution method', 'Sparsity pattern'
       else if ( ts_method == TS_SPARSITY_TRI ) then
          write(*,10)'Solution method', 'Sparsity pattern with tri-solver'
          if ( opt_TriMat_method == 0 ) then
             chars = 'speed'
          else if  ( opt_TriMat_method == 1 ) then
             chars = 'memory'
          end if
          write(*,10)'Tri-Mat create algorithm', trim(chars)
#ifdef MUMPS
       else if ( ts_method == TS_SPARSITY_MUMPS ) then
          write(*,10)'Solution method', 'Sparsity pattern + MUMPS library'
#endif
       end if
       write(*,8) 'TranSIESTA SCF cycle mixing weight',ts_wmix
       write(*,1) 'Start TS-SCF cycle immediately', ImmediateTSmode
       if ( IsVolt ) then
          write(*,6) 'Voltage', Volt/eV,'Volts'
          if ( VoltageInC ) then
             write(*,11) 'Voltage drop across central region'
          else
             write(*,11) 'Voltage drop across entire cell'    
          end if

          chars = 'Non-equilibrium contour weight method'
          select case ( TS_W_METHOD )
          case ( TS_W_ORB_ORB )
             write(*,10) trim(chars),'orb-orb'
          case ( TS_W_CORRELATED + TS_W_TR_ATOM_ATOM )
             write(*,10) trim(chars),'Correlated Tr[atom]-Tr[atom]'
          case ( TS_W_TR_ATOM_ATOM )
             write(*,10) trim(chars),'Uncorrelated Tr[atom]-Tr[atom]'
          case ( TS_W_CORRELATED + TS_W_TR_ATOM_ORB )
             write(*,10) trim(chars),'Correlated Tr[atom]-orb'
          case ( TS_W_TR_ATOM_ORB )
             write(*,10) trim(chars),'Uncorrelated Tr[atom]-orb'
          case ( TS_W_CORRELATED + TS_W_SUM_ATOM_ATOM )
             write(*,10) trim(chars),'Correlated Sum[atom]-Sum[atom]'
          case ( TS_W_SUM_ATOM_ATOM )
             write(*,10) trim(chars),'Uncorrelated Sum[atom]-Sum[atom]'
          case ( TS_W_CORRELATED + TS_W_SUM_ATOM_ORB )
             write(*,10) trim(chars),'Correlated Sum[atom]-orb'
          case ( TS_W_SUM_ATOM_ORB )
             write(*,10) trim(chars),'Uncorrelated Sum[atom]-orb'
          case default
             ! This is an easy place for cathing mistakes
             call die('Error in code, this should never occur.')
          end select
          chars = 'Non-equilibrium contour weight k-method'
          select case ( TS_W_K_METHOD ) 
          case ( TS_W_K_CORRELATED )
             write(*,10) trim(chars),'Correlated k-points'
          case ( TS_W_K_HALF_CORRELATED )
             write(*,10) trim(chars),'Half-correlated'
          case ( TS_W_K_UNCORRELATED )
             write(*,10) trim(chars),'Uncorrelated k-points'
          end select
       else
          write(*,11) 'TranSIESTA no voltage applied'
       end if
       if ( .not. Calc_Forces ) then
          write(*,11) '*** TranSIESTA will NOT update forces ***'
       end if

       if ( TS_RHOCORR_METHOD == 0 ) then
          write(*,11)'Will not correct charge fluctuations'
       else if ( TS_RHOCORR_METHOD == TS_RHOCORR_BUFFER ) then ! Correct in buffer
          if ( 0 < na_Buf ) then
             write(*,10)'Charge fluctuation correction','buffer'
          else
             call die('Charge correction can not happen in buffer as no buffer &
                  &atoms exist.')
          end if
          write(*,8)'Charge correction factor',TS_RHOCORR_FACTOR
       end if
       write(*,10)'          >> Electrodes << '
       do i = 1 , size(Elecs)
          write(*,11) '>> '//trim(name(Elecs(i)))
          if ( Elecs(i)%out_of_core ) then
             write(*,10) '  GF file', trim(GFFile(Elecs(i)))
             write(*,10) '  GF title', trim(GFtitle(Elecs(i)))
             write(*,1)  '  Reuse existing GF-file', Elecs(i)%ReUseGF
          else
             write(*,11)  '  In-core GF'
          end if
          write(*,10) '  Electrode TSHS file', trim(HSFile(Elecs(i)))
          write(*,5)  '  # atoms used in electrode', Elecs(i)%na_used
          write(*,15) '  Electrode repetition [A1 x A2 x A3]', &
               Elecs(i)%RepA1,Elecs(i)%RepA2,Elecs(i)%RepA3
          if ( Elecs(i)%t_dir == 1 ) then
             chars = 'A1'
          else if ( Elecs(i)%t_dir == 2 ) then
             chars = 'A2'
          else if ( Elecs(i)%t_dir == 3 ) then
             chars = 'A3'
          end if
          j = Elecs(i)%idx_na
          write(*,20) '  Position in geometry', j, j + TotUsedAtoms(Elecs(i)) - 1
          write(*,10) '  Transport direction for electrode', trim(chars)
          if ( Elecs(i)%inf_dir == INF_POSITIVE ) then
          write(*,10) '  Semi-infinite direction for electrode', 'positive wrt. '//trim(chars)
          else
          write(*,10) '  Semi-infinite direction for electrode', 'negative wrt. '//trim(chars)
          end if
          write(*,7)  '  Chemical shift', Elecs(i)%mu%mu/eV,'eV'
          write(*,1)  '  Bulk values in electrode', Elecs(i)%Bulk
          write(*,1)  '  Update cross terms contact/electrode', Elecs(i)%DM_CrossTerms
          write(*,1)  '  Calc. valence band-bottom eigenvalue', Elecs(i)%BandBottom
          write(*,8)  '  Hamiltonian E-C Ef fractional shift', Elecs(i)%Ef_frac_CT
       end do

       ! Print the contour information
       call print_contour_options( 'TS' , IsVolt )
       
       write(*,11) repeat('*', 62)
       write(*,*)

       write(*,'(3a)') repeat('*',24),' Begin: TS CHECKS AND WARNINGS ',repeat('*',24)

       ! Calculate the number of optimal contour points
       i = mod(N_Eq_E(), Nodes) ! get remaining part of equilibrium contour
       if ( IONode .and. i /= 0 ) then
          i = Nodes - i
          write(*,'(a)')'Without loosing performance you can increase &
               &the equilibrium integration precision.'
          write(*,'(a,i0,a)')'You can add ',i,' more energy points in the &
               &equilibrium contours, for FREE!'
          if ( i/N_mu > 0 ) then
             write(*,'(a,i0,a)')'This is ',i/N_mu,' more energy points per chemical potential.'
          end if
       end if

       i = mod(N_nEq_E(), Nodes) ! get remaining part of equilibrium contour
       if ( IONode .and. i /= 0 ) then
          i = Nodes - i
          write(*,'(a)')'Without loosing performance you can increase &
               &the non-equilibrium integration precision.'
          write(*,'(a,i0,a)')'You can add ',i,' more energy points in the &
               &non-equilibrium contours, for FREE!'
       end if

       if ( .not. Calc_Forces ) then
          write(*,11) '***       TranSIESTA will NOT update forces       ***'
          write(*,11) '*** ALL FORCES AFTER TRANSIESTA HAS RUN ARE WRONG ***'
       end if

       if ( ts_method == TS_SPARSITY_MUMPS .and. IsVolt ) then
          call die('Currently the bias contour is not functioning &
               &for the MUMPS solver.')
       end if

       ! Check that the unitcell does not extend into the transport direction
       do i = 1 , 3
          if ( i == ts_tdir ) cycle
          if ( abs(dot(ucell(:,i),ucell(:,ts_tdir),3)) > 1e-7_dp ) then
             write(*,*) &
                  "ERROR: Unitcell has the electrode extend into the &
                  &transport direction."
             write(*,*) &
                  "Please change the geometry."
             call die("Electrodes extend into the transport direction. &
                  &Please change the geometry.")
          end if
       end do

       ! Check that the atoms are placed correctly in the unit-cell
       ! The Hartree potential correction will only be put correctly 
       ! when the atoms are sorted by z and starting from z == 0
       if ( IsVolt .and. .not. VoltageInC ) then
          tmp = minval(xa(ts_tdir,:)) / Ang
          ! below -.5 or above .5 Ang from the bottom of the unit-cell
          if ( tmp < -.5_dp .or. .5_dp < tmp ) then
             write(*,*) &
                  "ERROR: Atoms must be located within the primary unit-cell &
                  &and shifted to 0 when dealing with bias."
             write(*,'(a,g15.6,a)') &
                  "Please shift your system in the z-direction by: ", &
                  -tmp,' Ang'
             call die('System setup wrong, please see output')
          end if
       end if

       ! warn the user about suspicous work regarding the electrodes
       do i = 1 , N_Elec

          if ( .not. Elecs(i)%Bulk ) then
             write(*,'(a)') 'Electrode '//trim(Name(Elecs(i)))//' will &
                  &not use bulk Hamiltonian. &
                  &Be careful here.'
          end if

          if ( Elecs(i)%DM_CrossTerms ) then
             write(*,'(a)') 'Electrode '//trim(Name(Elecs(i)))//' will &
                  &update cross-terms with central region.'
          end if

          if ( .not. Elecs(i)%kcell_check ) then
             write(*,'(a)') 'Electrode '//trim(Name(Elecs(i)))//' will &
                  &not check the k-grid sampling vs. system k-grid &
                  &sampling. Please ensure appropriate sampling.'
          end if
          if ( Elecs(i)%Ef_frac_CT /= 0._dp ) then
             write(*,'(a)') 'Electrode '//trim(Name(Elecs(i)))//' will &
                  &shift coupling Hamiltonian with a shift in energy. &
                  &Be careful here.'
          end if

       end do

    end if

    if ( IONode ) then
       write(*,'(3a,/)') repeat('*',24), &
            ' End: TS CHECKS AND WARNINGS ',repeat('*',26)
    end if


    if ( IONode ) then
       write(*,'(/,a,/)') '### Transiesta information for FDF-file START ###'
    end if
    
    call print_mus_block( 'TS' , N_mu , mus)

    call print_contour_block( 'TS' , IsVolt )

    if ( IONode ) then
       write(*,'(/,a,/)') '### Transiesta information for FDF-file END ###'
    end if

    ! write out the contour
    call io_contour(IsVolt, mus, kT, slabel)

    ! Print out the electrode coordinates
    do i = 1 , N_Elec
       call print_elec(Elecs(i),na_u,xa)
    end do
    
1   format('ts_options: ',a,t53,'=',4x,l1)
5   format('ts_options: ',a,t53,'=',i5,a)
20  format('ts_options: ',a,t53,'= ',i0,' -- ',i0)
6   format('ts_options: ',a,t53,'=',f10.4,tr1,a)
7   format('ts_options: ',a,t53,'=',f12.6,tr1,a)
8   format('ts_options: ',a,t53,'=',f10.4)
10  format('ts_options: ',a,t53,'=',4x,a)
11  format('ts_options: ',a)
15  format('ts_options: ',a,t53,'= ',i0,' x ',i0,' x ',i0)
    
  end subroutine read_ts_options
  
end module m_ts_options
