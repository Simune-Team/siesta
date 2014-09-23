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
  logical :: TS_HS_save = .true.
  logical :: TS_DE_save = .false.
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

  ! Controls how the initial guess for the density matrix will be formed
  ! Either it can be 'diagon' == 0, or 'transiesta' == 1 where the 
  ! latter is reading in the electrode DM and EDM
  integer :: TS_scf_mode = 0
  integer :: DM_bulk = 0

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
  
  subroutine read_ts_options( wmix, kT, ucell, Nmove, na_u, xa, lasto)

    use alloc
    use files, only : slabel
    use fdf, only : fdf_get, fdf_deprecated, fdf_obsolete
    use fdf, only : leqi
    use parallel, only: IOnode, Nodes
    use units, only: eV, Ang, Kelvin
    use intrinsic_missing, only : VNORM

    use m_ts_cctype
    use m_ts_global_vars, only : TSmode, ts_istep
    use m_ts_io, only : ts_read_TSHS_opt

    use m_ts_contour
    use m_ts_contour_eq,  only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E
#ifdef TRANSIESTA_WEIGHT_DEBUG
    use m_ts_contour_eq,  only : Eq_E, ID2idx, c2weight_eq
    use m_ts_contour_neq, only : nEq_E, N_nEq_ID, c2weight_neq, ID2mu, nEq_ID
#endif

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
    integer, intent(in) :: Nmove, na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)

! *******************
! * LOCAL variables *
! *******************
    real(dp) :: tmp
    logical :: err
    character(len=200) :: c, chars
    integer :: i, j, idx, idx1, idx2
#ifdef TRANSIESTA_WEIGHT_DEBUG
    integer, allocatable :: ID_mu(:)
    real(dp), allocatable :: rnID(:), rn(:), rw(:)
    type(ts_c_idx) :: cE
    complex(dp) :: W, ZW
#endif

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
    TS_HS_save = fdf_get('TS.HS.Save',.true.)
    TS_DE_save = fdf_get('TS.DE.Save',.false.)
    onlyS      = fdf_get('TS.onlyS',.false.)

    ! Read in the transport direction
    chars = fdf_get('TS.TransportDirection','c')
    if ( leqi(chars,'a') .or. leqi(chars,'a1') ) then
       ts_tdir = 1
    else if ( leqi(chars,'b') .or. leqi(chars,'a2') ) then
       ts_tdir = 2
    else if ( leqi(chars,'c') .or. leqi(chars,'a3') ) then
       ts_tdir = 3
    else if ( leqi(chars,'none') ) then
       ts_tdir = 0
    else
       call die('Transport direction not in [a|b|c|A1|A2|A3|none]')
    end if

    if ( TS_HS_save .and. FixSpin ) then
       write(*,*) 'Fixed spin not possible with Transiesta!'
       write(*,*) 'Electrodes with fixed spin is not possible with Transiesta !'
       call die('Stopping code')
    end if

    if ( onlyS .or. .not. TSmode ) then
       if ( IONode ) then
          write(*,1) 'Save H and S matrices', TS_HS_save
          write(*,1) 'Save DM and EDM matrices', TS_DE_save
          write(*,1) 'Save S and quit (onlyS)', onlyS
          write(chars,'(a,i0)') 'A',ts_tdir
          write(*,10) 'Transport along unit-cell vector',trim(chars)
          write(*,11) repeat('*', 62)
          write(*,*)
       end if
       return
    end if

    ! Read in the mixing for the transiesta cycles
    ts_wmix = fdf_get('TS.MixingWeight',wmix)
    
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
    MUMPS_block = fdf_get('TS.MUMPS.BlockingFactor',112)
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
    chars = fdf_get('SCF.Initialize','diagon')
    if ( leqi(chars,'diagon') ) then
       TS_scf_mode = 0
    else if ( leqi(chars,'transiesta') ) then
       TS_scf_mode = 1
    end if

    ! Whether we should always set the DM to bulk
    ! values (by reading in from electrode DM)
    chars = fdf_get('TS.Elecs.DM.Bulk','none')
    DM_bulk = 0
    if ( leqi(chars,'init') ) then
       DM_bulk = 1
    !else if ( leqi(chars,'scf') ) then
    !   DM_bulk = 2
    end if

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
    else if ( leqi(chars,'fermi') ) then
       TS_RHOCORR_METHOD = TS_RHOCORR_FERMI
    end if
    TS_RHOCORR_FERMI_TOLERANCE = &
         fdf_get('TS.ChargeCorrection.Fermi.Tolerance',0.01_dp)
    ! Factor for charge-correction
    TS_RHOCORR_FACTOR = fdf_get('TS.ChargeCorrection.Factor',0.75_dp)
    if ( TS_RHOCORR_METHOD == TS_RHOCORR_BUFFER ) then
       if ( 1.0_dp < TS_RHOCORR_FACTOR ) then
          call die("Charge correction factor must be in the range [0;1]")
       endif
    end if
    if ( TS_RHOCORR_FACTOR < 0.0_dp ) then
       call die("Charge correction factor must be larger than 0")
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
       Elecs(:)%DM_update = 2
    else
       ! default to not update the cross-terms
       c = fdf_get('TS.Elecs.DM.Update','none')
       if ( leqi(c,'none') ) then
          Elecs(:)%DM_update = 0
       else if ( leqi(c,'cross-terms') ) then
          Elecs(:)%DM_update = 1
       else if ( leqi(c,'all') ) then
          Elecs(:)%DM_update = 2
       else
          call die('TS.Elecs.DM.Update [none,cross-terms,all]: &
               &unrecognized option: '//trim(c))
       end if
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
       if ( Elecs(i)%idx_a < 0 ) &
            Elecs(i)%idx_a = na_u + Elecs(i)%idx_a + 1
       if ( Elecs(i)%idx_a < 1 .or. &
            na_u < Elecs(i)%idx_a ) &
            call die("Electrode position does not exist")
       Elecs(i)%idx_o = lasto(Elecs(i)%idx_a-1)+1

    end do

    ! Check that we can actually start directly in transiesta
    if ( TS_scf_mode == 1 ) then ! TS-start
       if ( .not. all(Elecs(:)%DM_update >= 1) ) then
          call die('Requesting immediate start, yet we do not update &
               &cross-terms.')
       end if
    end if

    ! Check that the current transport direction is "aligned"
    ! TODO when we can deal with arbitrary electrodes this should be altered
    if ( N_Elec == 2 ) then
       ! For the moment we require that 2 electrodes have 
       ! a well defined transport direction
       if ( all(Elecs(1)%t_dir == Elecs(:)%t_dir) ) then
          ts_tdir = Elecs(1)%t_dir
       else
          call die('Using 2 electrodes requires a well-defined &
               &transport direction. Only use one direction for &
               the semi-infinite leads.')
       end if
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
       idx1 = Elecs(i)%idx_a
       idx2 = idx1 + TotUsedAtoms(Elecs(i)) - 1
       ! we need to check every electrode,
       ! specifically because if one of the electrodes is fully located
       ! inside the other and we check the "small" one 
       do j = 1 , N_Elec
          if ( i == j ) cycle
          idx = Elecs(j)%idx_a
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
             idx1 = Elecs(j)%idx_a
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
    if ( Nmove > 0 .and. .not. all(Elecs(:)%DM_update > 0) ) then
       call die('transiesta relaxation is only allowed if you also &
            &update the cross terms, please set: TS.Elecs.DM.Update cross-terms')
    end if
    if ( Nmove > 0 .and. .not. Calc_Forces ) then
       call die('transiesta relaxation is based on calculating the forces, &
            &you cannot relax your system and not require the calculation of forces')
    end if

    ! Update the weight function
    chars = fdf_get('TS.Weight.k.Method','correlated')
    if ( leqi(chars,'correlated') ) then
       TS_W_K_METHOD = TS_W_K_CORRELATED
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
               &must be [[un]correlated+][orb-orb|tr-atom-atom|sum-atom-atom|mean]')
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
    else if ( leqi(chars,'mean') ) then
       TS_W_METHOD = TS_W_MEAN
    else
       call die('Unrecognized option for TS.Weight.Method &
            &must be [[un]correlated+|][orb-orb|tr-atom-[atom|orb]|sum-atom-[atom|orb]|mean]')
    end if

    ! read in contour options
    if ( TSmode ) then
       call read_contour_options( N_Elec, Elecs, N_mu, mus, kT, IsVolt, Volt )
    end if

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
       write(*,1) 'Save H and S matrices', TS_HS_save
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
       if ( ts_tdir < 1 ) then
          write(*,11) 'Transport individually selected for electrodes'
       else
          write(chars,'(a,i0)') 'A',ts_tdir
          write(*,10) 'Transport along unit-cell vector',trim(chars)
       end if
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

       select case ( TS_scf_mode )
       case ( 0 )
          write(*,10) 'Initialize DM by','diagon'
       case ( 1 )
          write(*,10) 'Initialize DM by','transiesta'
       end select
       select case ( DM_bulk ) 
       case ( 0 ) 
          write(*,11) 'DM for electrodes will not be used'
       case ( 1 )
          write(*,11) 'DM for electrodes will be initialized to bulk'
       end select
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
          case ( TS_W_MEAN )
             write(*,10) trim(chars),'Algebraic mean'
          case default
             ! This is an easy place for cathing mistakes
             call die('Error in code, weighting method unrecognized.')
          end select
          chars = 'Non-equilibrium contour weight k-method'
          select case ( TS_W_K_METHOD ) 
          case ( TS_W_K_CORRELATED )
             write(*,10) trim(chars),'Correlated k-points'
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
       else if ( TS_RHOCORR_METHOD == TS_RHOCORR_FERMI ) then ! Correct fermi-lever
          write(*,10)'Charge correction','Fermi-level'
          write(*,8)'Charge correction tolerance',TS_RHOCORR_FERMI_TOLERANCE
          write(*,8)'Charge correction factor',TS_RHOCORR_FACTOR

       end if
       write(*,10)'          >> Electrodes << '
       do i = 1 , size(Elecs)
          write(*,11) '>> '//trim(name(Elecs(i)))
          if ( Elecs(i)%out_of_core ) then
             write(*,10) '  GF file', trim(Elecs(i)%GFfile)
             write(*,1)  '  Reuse existing GF-file', Elecs(i)%ReUseGF
          else
             write(*,11)  '  In-core GF'
          end if
          write(*,10) '  Electrode TSHS file', trim(Elecs(i)%HSfile)
          write(*,5)  '  # atoms used in electrode', Elecs(i)%na_used
          write(*,15) '  Electrode repetition [A1 x A2 x A3]', Elecs(i)%Rep(:)
          if ( Elecs(i)%t_dir == 1 ) then
             chars = 'A1'
          else if ( Elecs(i)%t_dir == 2 ) then
             chars = 'A2'
          else if ( Elecs(i)%t_dir == 3 ) then
             chars = 'A3'
          end if
          j = Elecs(i)%idx_a
          write(*,20) '  Position in geometry', j, j + TotUsedAtoms(Elecs(i)) - 1
          write(*,10) '  Transport direction for electrode', trim(chars)
          if ( Elecs(i)%inf_dir == INF_POSITIVE ) then
          write(*,10) '  Semi-infinite direction for electrode', 'positive wrt. '//trim(chars)
          else
          write(*,10) '  Semi-infinite direction for electrode', 'negative wrt. '//trim(chars)
          end if
          write(*,7)  '  Chemical shift', Elecs(i)%mu%mu/eV,'eV'
          write(*,1)  '  Bulk values in electrode', Elecs(i)%Bulk
          if ( product(Elecs(i)%Rep) > 1 ) then
             write(*,1)  '  Pre-expansion to reduce computation', Elecs(i)%pre_expand
          end if
          if ( Elecs(i)%DM_update == 0 ) then
             write(*,11)  '  Cross-terms is not updated'
          else if ( Elecs(i)%DM_update == 1 ) then
             write(*,11)  '  Cross-terms is updated'
          else if ( Elecs(i)%DM_update == 2 ) then
             write(*,11)  '  Cross-terms and electrode region is updated'
          end if
          write(*,1)  '  Calc. valence band-bottom eigenvalue', Elecs(i)%BandBottom
          if ( IsVolt ) &
               write(*,8)  '  Hamiltonian E-C Ef fractional shift', Elecs(i)%Ef_frac_CT
          if ( .not. Elecs(i)%kcell_check ) then
             write(*,11)  '  Will NOT check the kgrid-cell! Ensure sampling!'
          end if
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

       if ( ts_tdir == 0 ) then
          write(*,11) '*** TranSIESTA transport direction is arbitrary  ***'
       else if ( ts_tdir < 0 ) then
          write(*,11) '*** TranSIESTA transport direction is individual ***'
       end if

       ! Check that the unitcell does not extend into the transport direction
       do i = 1 , 3
          if ( i == ts_tdir .or. ts_tdir <= 0 ) cycle
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

       ! If the user has requested to initialize using transiesta
       ! and the user does not utilize the bulk DM, they should be
       ! warned
       if ( TS_scf_mode == 1 .and. DM_bulk == 0 ) then
          write(*,'(a)') 'You are not initializing the electrode DM/EDM. &
               &This may result in very wrong electrostatic potentials close to &
               &the electrode/device boundary region.'
       end if
          

       ! warn the user about suspicous work regarding the electrodes
       do i = 1 , N_Elec

          if ( .not. Elecs(i)%Bulk ) then
             write(*,'(a)') 'Electrode '//trim(Name(Elecs(i)))//' will &
                  &not use bulk Hamiltonian. &
                  &Be careful here.'
          end if

          if ( Elecs(i)%DM_update == 0 ) then
             write(*,'(a)') 'Electrode '//trim(Name(Elecs(i)))//' will &
                  &not update cross-terms or local region.'
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

          ! if any buffer atoms exist, we should suggest to the user
          ! to use TS.Elec.<elec> [DM-update cross-terms|all]
          ! in case any buffer atoms are too close
          err = .false.
          do j = 1 , na_u
             if ( atom_type(j) /= TYP_BUFFER ) cycle
             do idx = 0 , TotUsedAtoms(Elecs(i)) - 1
                ! Proximity of 12 Bohr enables this check
                err = vnorm(xa(:,Elecs(i)%idx_a+idx)-xa(:,j)) < 12._dp
                if ( err ) exit
             end do
             if ( err ) exit
          end do
          if ( err .and. Elecs(i)%DM_update == 0 ) then
             ! some buffer atoms are close to this electrode
             ! Advice to use dm_update
             write(*,'(a,/,a)') 'Electrode '//trim(Name(Elecs(i)))//' is &
                  &likely terminated by buffer atoms. It is HIGHLY recommended to add this:', &
                  '  TS.Elec.'//trim(Name(Elecs(i)))//' DM-update [cross-terms|all]'
          end if

          ! In case DM_bulk is requested we assert that the file exists
          inquire(file=Elecs(i)%DEfile,exist=err)
          err = .not. err
          if ( DM_bulk == 1 .and. err ) then
             write(*,'(a,/,a)') 'Electrode '//trim(Name(Elecs(i)))//' TSDE &
                  &file cannot be located in: '//trim(Elecs(i)%DEfile)//'.', &
                  '  Please add TS.DE.Save T to the electrode calculation or &
                  &specify the exact file position using ''TSDE-file'' in the&
                  & Elec block.'
          end if
          
       end do

       if ( N_Elec /= 2 .and. any(Elecs(:)%DM_update == 0) ) then
          write(*,'(a,/,a)') 'Consider updating more elements when doing &
               &N-electrode calculations. The charge conservation typically &
               &increases.','  TS.Elecs.DM.Update [cross-terms|all]'
       end if

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

#ifdef TRANSIESTA_WEIGHT_DEBUG
    if ( IONode ) then
       allocate(ID_mu(N_nEq_ID))
       allocate(rnID(N_nEq_ID),rw(N_mu),rn(N_mu))
       do i = 1 , N_nEq_ID
          ID_mu(i) = ID2mu(i)
          rnID(i) = i
       end do
       write(*,'(a)') 'Equilibrium:'
       tmp = .5_dp / 3.14159265358979323846_dp
       i = 1
       cE = Eq_E(i,step=1) ! we read them backwards
       do while ( cE%exist ) 
          
          do j = 1 , N_mu
             if ( cE%fake ) cycle
             call ID2idx(cE,mus(j)%ID,idx)
             if ( idx < 1 ) cycle
             call c2weight_eq(cE,idx, tmp, W ,ZW)

             write(*,'(i2,tr1,a10,2(tr1,i2),4(tr1,f10.5))') &
                  i,trim(mus(j)%name),mus(j)%ID,idx,W,ZW / eV
          end do
          i = i + 1
          cE = Eq_E(i,step=1)
       end do

       write(*,'(a)') 'Non-equilibrium:'
       i = 1
       cE = nEq_E(i,step=1) ! we read them backwards
       do while ( cE%exist ) 
          
          do j = 1 , N_Elec
             if ( cE%fake ) cycle
             if ( .not. has_cE(cE,iEl=j) ) cycle

             do idx1 = 1 , N_nEq_ID
                if ( .not. has_cE(cE,iEl=j,ineq=idx1) ) cycle
                
                call c2weight_neq(cE,kT,j,idx1, tmp,W,idx,ZW)

                write(*,'(i2,tr1,a10,2(tr1,i2),4(tr1,f10.5))') &
                     i,trim(Elecs(j)%name),mus(idx)%ID,idx1,W,ZW / eV
             end do
          end do

          i = i + 1
          cE = nEq_E(i,step=1)
       end do
       write(*,'(a)') 'DM_neq: '
       do i = 1 , N_nEq_ID
          write(*,'(2(a10,tr1),f10.5)') &
               trim(Elecs(nEq_ID(i)%iEl)%name),trim(mus(nEq_ID(i)%imu)%name),rnID(i)
       end do
       call calc_neq_weight(N_Elec,N_mu,N_nEq_ID,ID_mu,rnID,rn,rw)
       write(*,'(a)') 'Contrib and weights: '
       do i = 1 , N_mu
          write(*,'(a10,2(tr1,f10.5))') trim(mus(i)%name),rn(i),rw(i)
       end do
       !call die('Stopping on request! Debugging WEIGHTS!!!')
    end if
#endif


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
