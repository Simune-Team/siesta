!==========================================================================*
!                                                                          *
!  TRANSIESTA MODULE m_ts_options : Declaration of the variables           *
!  involved in a TRANSIESTA calculation.                                   *
!                                                                          *
!  Written by F.D.Novaes, May'07                                           *
!  onlyS option added by M.Paulsson May'09                                 *
!==========================================================================*
!  Contains the Subroutines:                                               *
!                                                                          *
!  1) read_ts_options : Reads the optional parameters from the fdf file    *
!                                                                          *
!==========================================================================* 


module m_ts_options

! SIESTA Modules used
  USE precision, only : dp
  USE siesta_options, only : FixSpin, isolve, SOLVE_TRANSI
  USE sys, only : die
  USE m_ts_electype
  use m_ts_chem_pot
  use m_ts_tdir
  implicit none
  PUBLIC
  SAVE

!=========================================================================*
!  Arguments read from input file using the fdf package                    *
!--------------------------------------------------------------------------*
  
logical  :: SaveTSHS = .true.     ! Saves the Hamiltonian and Overlap matrices if the 
                         ! the option TS.SaveHS is specified in the input file
logical  :: onlyS = .false. ! Option to only save overlap matrix
logical  :: IsVolt = .false.      ! Logical for dabs(VoltFDF) > 0.0001d*eV
real(dp) :: Volt = 0._dp         ! Bias applied, Internally Volt=voltfdf/eV (eV). 
integer  :: na_BufL = 0     ! Number of Left Buffer Atoms
integer  :: na_BufR = 0     ! Number of Right Buffer Atoms
integer  :: no_BufL = 0     ! Number of Left Buffer orbitals
integer  :: no_BufR = 0     ! Number of Right Buffer orbitals
! Electrodes and different chemical potentials
integer :: N_Elec = 0
type(Elec), allocatable, target :: Elecs(:)
integer :: N_mu = 0
type(ts_mu), allocatable, target :: mus(:)
logical :: ReUseGF = .false.         ! Calculate the electrodes GF
logical :: ImmediateTSmode = .false. ! will determine to immediately start the transiesta
                                     ! SCF. This is useful when you already have a converged
                                     ! siesta DM

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

real(dp) :: Elecs_xa_EPS = 1.e-4_dp

! The mixing weight in the transiesta cycles...
real(dp) :: ts_wmix ! = wmix

! The user can request to analyze the system, returning information about the 
! tri-diagonalization partition and the contour
logical :: TS_Analyze = .false.
integer :: TS_bandwidth_algo = 0

! Flag to control TranSIESTA
logical :: TSmode = .false.

CONTAINS

! *********************************************************************
! Subroutine to read the data for the TRANSIESTA program
!
!     It uses the FDF (Flexible Data Format) package
!     of J.M.Soler and A.Garcia
!
! Writen by F.D.Novaes May'07
! Rewritten by Nick Papior Andersen, 2013
!
! **************************** OUTPUT *********************************
  
  subroutine read_ts_options( wmix, kT, ucell, na_u, xa, lasto)

! SIESTA Modules Used
    use alloc
    use files, only  : slabel
    use fdf, only : fdf_get, fdf_deprecated, fdf_obsolete
    use fdf, only : leqi
    use parallel, only: IOnode, Nodes, operator(.parcount.)
    use units, only: eV, Ang, Kelvin
    use m_ts_cctype
    use m_ts_global_vars, only : ts_istep, TSinit
    use m_ts_io, only : ts_read_TSHS_opt

    use m_ts_contour
    use m_ts_io_contour
    use m_ts_method
    use m_ts_weight
    use m_ts_charge
    use m_ts_tdir
    
    use m_monitor
    use m_bandwidth

    implicit none
    
    real(dp), intent(in) :: wmix, kT
    real(dp),intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
! Internal Variables
    real(dp) :: tmp
    logical :: err
    character(len=200) :: c, chars
    integer :: i, j, idx, idx1, idx2

    ! External routines
    real(dp) :: dot
    external :: dot

    if (isolve.eq.SOLVE_TRANSI) then
       TSmode = .true.
       ! If in TSmode default to initalization
       ! In case of 'DM.UseSaveDM TRUE' TSinit will be set accordingly
       TSinit = .true.
    endif

    if (IOnode) then
       write(*,*)
       write(*,11) repeat('*', 62)
    end if

    !Set ts_istep default
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
    else
       call die('Unrecognized Transiesta solution method: '//trim(chars))
    end if


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
       if ( leqi(chars,'cuthill-mckee') ) then
          TS_bandwidth_algo = TS_bandwidth_algo + BW_CUTHILL_MCKEE
       else if ( leqi(chars,'cuthill-mckee-z-priority') .or. &
            leqi(chars,'cuthill-mckee-z') ) then
          TS_bandwidth_algo = TS_bandwidth_algo + BW_CUTHILL_MCKEE_Z_PRIORITY
       else if ( leqi(chars,'papior') ) then
          TS_bandwidth_algo = TS_bandwidth_algo + BW_PAPIOR
       else
          call die('Unrecognized option for Bandwidth algorithm: '//trim(chars))
       end if
    end if
    
    ! Update the weight function
    chars = fdf_get('TS.Weight.NonEquilibrium','correlated')
    if ( leqi(chars,'correlated') ) then
       TS_W_METHOD = TS_W_CORRELATED
    else if ( leqi(chars,'uncorrelated') ) then
       TS_W_METHOD = TS_W_UNCORRELATED
       call die('Not currently functioning')
    else if ( leqi(chars,'k-uncorrelated') ) then
       TS_W_METHOD = TS_W_K_UNCORRELATED
    else
       call die('Could not determine flag TS.Weight.NonEquilibrium, please &
            &see manual.')
    end if

    ! Figure out the number of orbitals on the buffer atoms
    na_BufL = fdf_get('TS.BufferAtomsLeft',0)
    no_BufL = 0
    do i = 1 , na_BufL
       no_BufL = no_BufL + lasto(i) - lasto(i-1)
    end do

    na_BufR = fdf_get('TS.BufferAtomsRight',0)
    no_BufR = 0
    do i = na_u - na_BufR + 1 , na_u
       no_BufR = no_BufR + lasto(i) - lasto(i-1)
    end do

    ! check that it is correctly setup
    if ( na_BufL < 0 .or. na_BufR < 0 ) then
       call die("Buffer atoms must be 0 or a positive integer.")
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
    
    Calc_Forces = fdf_get('TS.Forces.Calc',.true.)

    call fdf_deprecated('TS.CalcGF','TS.ReUseGF')
    ReUseGF = fdf_get('TS.ReUseGF',.false.)

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
    if ( N_mu < 1 ) &
         call die('You need at least one chemical potential')
    do i = 1 , N_mu
       ! Default things that could be of importance
       if ( .not. fdf_mu('TS',slabel,mus(i)) ) then
          call die('Could not find chemical potential: '//trim(name(mus(i))))
       end if
       ! Attach the ID
       mus(i)%ID = i
    end do

    ! To determine the same coordinate nature of the electrodes
    Elecs_xa_EPS= fdf_get('TS.Elecs.Coord.Eps',1.e-4_dp,'Bohr')

    ! detect how many electrodes we have
    N_Elec = fdf_nElec('TS',Elecs)
    if ( N_Elec < 2 ) then
       call die('Please see the manual for how to construct an &
            &example electrode configuration')
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
    Elecs(:)%BandBottom = fdf_get('TS.Elecs.Calc.BandBottom', .false.)

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
    if ( all(Elecs(1)%t_dir == Elecs(:)%t_dir) ) then
       ts_tdir = Elecs(1)%t_dir
    end if

    ! The sign can not be chosen from this (several mu, where to define it)
    Volt = maxval(mus(:)%mu) - minval(mus(:)%mu)
    call fdf_obsolete('TS.Voltage')
write(*,*) 'TODO the bias is not determined correctly by the direction, see m_ts_voltage, say if we change sign'
    !Volt     = fdf_get('TS.Voltage',0._dp,'Ry') 
    ! Voltage situation is above 0.01 mV
    IsVolt = dabs(Volt) > 0.00001_dp*eV
    if ( .not. IsVolt ) then
       Volt = 0._dp
       mus(:)%mu = 0._dp

       ! We must make all electrodes point to the first chemical potential 
       ! and discard the rest

       ! Save che first chemical potential name
       c = mus(1)%name
       do i = 1 , N_mu
          deallocate(mus(i)%Eq_seg)
       end do
       deallocate(mus)
       ! create the first chemical potential again
       N_mu = 1
       allocate(mus(1))
       mus(1)%name = trim(c)
       if ( .not. fdf_mu('TS',slabel,mus(1)) ) then
          call die('Could not find chemical potential: '//trim(name(mus(1))))
       end if
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
          if ( abs(mus(i)%mu - mus(j)%mu) > 1.e-3_dp*eV ) then
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
       if ( idx1 <= na_BufL ) then
          write(*,*) 'Electrode: '//trim(Name(Elecs(i)))
          write(*,'(a,i0,a,i0)') 'Positions: ',idx1,' -- ',idx2 
          write(*,'(a,i0)') 'Buffer atoms stops at: ',na_BufL
          call die('Left buffer atoms overlap an electrode')
       else if ( na_u - na_BufR < idx2 ) then
          write(*,*) 'Electrode: '//trim(Name(Elecs(i)))
          write(*,'(a,i0,a,i0)') 'Positions: ',idx1,' -- ',idx2 
          write(*,'(a,i0)') 'Buffer atoms starts at: ',na_u - na_BufR + 1
          call die('Right buffer atoms overlap an electrode')
       end if
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
          idx = Elecs(j)%idx_na + TotUsedAtoms(Elecs(j)) - 1
          if ( (idx <= idx1 .and. &
               idx1 < idx + TotUsedAtoms(Elecs(j))) ) then
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
    if (  .5_dp * abs(Volt) < maxval(mus(:)%mu) .or. &
         -.5_dp * abs(Volt) > minval(mus(:)%mu) ) then
       if ( IONode ) then
          write(*,'(a)') 'Chemical potentials [eV]:'
          do i = 1 , N_Elec
             write(*,'(a,f10.5,a)') trim(Name(Elecs(i)))//' at ',Elecs(i)%mu%mu/eV,' eV'
          end do
          write(*,'(a,f10.5,a)') '-V/2'//' at ',-Volt*.5_dp/eV,' eV'
          write(*,'(a,f10.5,a)') ' V/2'//' at ',Volt*.5_dp/eV,' eV'
       end if
       call die('Chemical potentials are not in range [-V/2 ; V/2]')
    end if

    ! WILL WORK EVENTUALLY
    write(*,*) 'TODO SEVERAL ELECTRODES POTENTIAL DROP!'
    !if ( size(Elecs) > 2 ) call die('currently does not work')

    ! read in contour options
    if ( TSmode ) then
       call read_contour_options( N_Elec, Elecs, N_mu, mus, kT, IsVolt, Volt )
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
          case ( BW_CUTHILL_MCKEE_Z_PRIORITY ) 
             chars = trim(chars)//'Cuthill-Mckee with z-priority'
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
          chars = 'Non-equilibrium contour weights'
          if  ( TS_W_METHOD == TS_W_CORRELATED ) then
             write(*,10) trim(chars),'Real-space'
          else if ( TS_W_METHOD == TS_W_UNCORRELATED ) then
             write(*,10) trim(chars),'Uncorrelated k-points in real-space'
          else if ( TS_W_METHOD == TS_W_K_UNCORRELATED ) then
             write(*,10) trim(chars),'Uncorrelated k-points in k-space'
          end if
       else
          write(*,11) 'TranSIESTA no voltage applied'
       end if
       if ( .not. Calc_Forces ) then
          write(*,11) '*** TranSIESTA will NOT update forces ***'
       end if
       write(*,5) 'Left buffer atoms', na_BufL
       write(*,5) 'Right buffer atoms', na_BufR


       if ( TS_RHOCORR_METHOD == 0 ) then
          write(*,11)'Will not correct charge fluctuations'
       else if ( TS_RHOCORR_METHOD == TS_RHOCORR_BUFFER ) then ! Correct in buffer
          if ( 0 < na_BufL .or. 0 < na_BufR ) then
             write(*,10)'Charge fluctuation correction','buffer'
          else
             call die('Charge correction can not happen in buffer as no buffer &
                  &atoms exist.')
          end if
          write(*,8)'Charge correction factor',TS_RHOCORR_FACTOR
       end if
       write(*,1) 'Re-use GF file if exists', ReUseGF
       write(*,10)'          >> Electrodes << '
       do i = 1 , size(Elecs)
          write(*,11) '>> '//trim(name(Elecs(i)))
          write(*,10) '  GF file', trim(GFFile(Elecs(i)))
          write(*,10) '  GF title', trim(GFtitle(Elecs(i)))
          write(*,10) '  Electrode TSHS file', trim(HSFile(Elecs(i)))
          write(*,5)  '  # atoms used in electrode', UsedAtoms(Elecs(i))
          write(*,15) '  Electrode repetition [A1 x A2 x A3]', &
               RepA1(Elecs(i)),RepA2(Elecs(i)),RepA3(Elecs(i))
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
       end do

       ! Print the contour information
       call print_contour_options( 'TS' , IsVolt )
       
       write(*,11) repeat('*', 62)
       write(*,*)

       write(*,'(3a)') repeat('*',24),' Begin: TS CHECKS AND WARNINGS ',repeat('*',24)

       if ( .not. Calc_Forces ) then
          write(*,11) '*** TranSIESTA will NOT update forces ***'
       end if


       ! Check that the unitcell does not extend into the transport direction
       do i = 1 , 3
          if ( i == ts_tdir ) cycle
          if ( abs(dot(ucell(:,i),ucell(:,ts_tdir),3)) > 1e-7 ) then
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
write(*,*)'Check the bias placement'
       if ( IsVolt .and. .false. ) then
          tmp = minval(xa(3,na_BufL+1:)) / Ang
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

    end if

    if ( IONode ) then
       write(*,'(3a,/)') repeat('*',24), &
            ' End: TS CHECKS AND WARNINGS ',repeat('*',26)
    end if

    call print_contour_block( 'TS' , IsVolt )

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
