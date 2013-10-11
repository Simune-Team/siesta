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
  USE siesta_options, only : fixspin, isolve, SOLVE_TRANSI
  USE sys, only : die
  USE m_ts_electype
  use m_ts_tdir
  implicit none
  PUBLIC
  SAVE

!=========================================================================*
!  Arguments read from input file using the fdf package                    *
!--------------------------------------------------------------------------*
  
logical  :: savetshs     ! Saves the Hamiltonian and Overlap matrices if the 
                         ! the option TS.SaveHS is specified in the input file
logical  :: onlyS        ! Option to only save overlap matrix
logical  :: USEBULK      ! Use Bulk Hamiltonian in Electrodes
  ! Logical variable that describes the solution method on
  ! LEFT-RIGHT-EQUILIBRIUM contour points.
  ! This is realized by the fact that for:
  !    UseBulk .and. UpdateDMCR
  ! the needed part of the GF is only the C...C regions:
  !
  !  -------------------------------------
  !  | L...L | L...C   0     0   |   0   |
  !  | C...L | C...C C...C C...C |   0   |
  !  |   0   | C...C C...C C...C |   0   |
  !  |   0   | C...C C...C C...C | C...R |
  !  |   0   |   0     0   R...C | R...R |
  !  -------------------------------------
  ! 
  ! This means we can solve the following instead:
  ! G_F^{-1} G_F I_P = I \times I_P,
  ! where I_P:
  !  ---------------------
  !  |   0     0     0   |
  !  |   1     0     0   |
  !  |   0     1     0   |
  !  |   0     0     1   |
  !  |   0     0     0   |
  !  ---------------------
  ! Note, that this can ONLY be used in EQUI contour points.
  ! In principle we can obtain the EXACT size of the problem
  ! For very large electrodes. This could come in handy.
logical  :: UpdateDMCR   ! Update DM values of ONLY Central Region
logical  :: UseVFix      ! Call the routine TSVHFix 
logical  :: IsVolt       ! Logical for dabs(VoltFDF) > 0.001d/eV
real(dp) :: VoltFDF      ! Bias applied, Internally Volt=voltfdf/eV. 
real(dp) :: VoltL        ! Bias on the left electrode   (  .5 * VoltFDF )
real(dp) :: VoltR        ! Bias on the right electrode  ( -.5 * VoltFDF )
integer  :: na_BufL      ! Number of Left Buffer Atoms
integer  :: na_BufR      ! Number of Right Buffer Atoms
integer  :: no_BufL      ! Number of Left Buffer orbitals
integer  :: no_BufR      ! Number of Right Buffer orbitals
type(Elec), allocatable :: Elecs(:)
logical :: ElecValenceBandBot ! Calculate Electrode valence band bottom when creating electrode GF
logical :: ReUseGF        ! Calculate the electrodes GF
logical :: ImmediateTSmode=.false. ! will determine to immediately start the transiesta
                           ! SCF. This is useful when you already have a converged
                           ! siesta DM

! If the energy-contour is not perfectly divisable by the number of nodes then adjust
integer :: opt_TriMat_method = 0 ! Optimization method for determining the best tri-diagonal matrix split
! 0  == We optimize for speed
! 1  == We optimize for memory

logical :: VoltageInC ! Determines whether the voltage-drop should be located in the constriction
                      ! I.e. if the electrode starts at 10 Ang and the central region ends at 20 Ang
                      ! then the voltage drop will only take place between 10.125 Ang and 19.875 Ang

real(dp) :: Elec_xa_EPS

! The mixing weight in the transiesta cycles...
real(dp) :: ts_wmix

! The user can request to analyze the system, returning information about the 
! tri-diagonalization partition and the contour
logical, save :: TS_Analyze = .false.
integer, save :: TS_bandwidth_algo = 0

! If the user request to monitor the Density matrix update elements
integer,          save :: N_mon = 0
integer, pointer, save :: monitor_list(:,:) => null()
integer, pointer, save :: iu_MON(:,:) => null()

!==========================================================================*
!==========================================================================*
!  Default Values for arguments read from input file                       *
!--------------------------------------------------------------------------*

logical, parameter :: savetshs_def = .true.
logical, parameter :: onlyS_def = .false.
logical, parameter :: tsdme_def = .true.
logical, parameter :: UseBulk_def = .true.
logical, parameter :: UpdateDMCR_def = .true.
logical, parameter :: UseVFix_def = .true.
real(dp), parameter :: voltfdf_def = 0._dp   ! in Ry
integer, parameter :: na_BufL_def = 0
integer, parameter :: na_BufR_def = 0
integer, parameter :: NRepA_def = 1
integer, parameter :: NUsedAtoms_def = -1
character(20), parameter :: smethod_def = 'gaussfermi'
character(33), parameter :: GFTitle_def = 'Generated GF file'
character(33), parameter :: HSFile_def = 'NOT REQUESTED'
logical, parameter :: ElecValenceBandBot_def = .false.
logical, parameter :: ReUseGF_def = .false.

logical, save :: TSmode = .false.

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
  
  subroutine read_ts_options( wmix, kT_in, ucell, na_u, xa, lasto)

! SIESTA Modules Used
    use files, only  : slabel
    use fdf, only : fdf_get, fdf_deprecated, fdf_obsolete
    use fdf, only : leqi
    use parallel, only: IOnode, Nodes, operator(.parcount.)
    use units, only: eV, Ang
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
    
    real(dp), intent(in) :: wmix, kT_in
    real(dp),intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
! Internal Variables
    real(dp) :: tmp
    character(len=70) :: chars
    integer :: tmp_G_NF
    integer :: i, j, idx, idx1, idx2
    type(Elec) :: tmpElec


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

    ! Read in the mixing for the transiesta cycles
    ts_wmix = fdf_get('TS.MixingWeight',wmix)
    
    ! Read in the transport direction
    ts_tdir = fdf_get('TS.TransportDirection',3)
    if ( ts_tdir < 1 .or. 3 < ts_tdir ) then
       call die('Transport direction not in [1-3]')
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

    ! Reading TS Options from fdf ...
    saveTSHS    = fdf_get('TS.SaveHS',savetshs_def)
    onlyS       = fdf_get('TS.onlyS',onlyS_def)

    VoltFDF     = fdf_get('TS.Voltage',voltfdf_def,'Ry') 
    ! Voltage situation is above 1 meV (probably too low...)
    IsVolt = dabs(VoltFDF) > 0.001_dp*eV
    if ( .not. IsVolt ) VoltFDF = 0._dp

    ! Set up the fermi shifts for the left and right electrodes
    VoltL =  0.5_dp * VoltFDF
    VoltR = -0.5_dp * VoltFDF

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

    ! currently this does not work
    !ImmediateTSmode = fdf_get('TS.SCFImmediate',.false.)

    UseBulk    = fdf_get('TS.UseBulkInElectrodes',UseBulk_def)
    UpdateDMCR = fdf_get('TS.UpdateDMCROnly',UpdateDMCR_def)
    if ( .not. UseBulk ) then
       ! If we use bulk, the algorithms look up in the UpdateDMCR
       ! to check whether the off-diagonal elements of
       ! the Green's function is needed
       ! Hence if we update everything, everything needs to
       ! be calculated: UpdateDMCR == .false.
       UpdateDMCR = .false.
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
    
    ! Update the weight function
    chars = fdf_get('TS.Weight.NonEquilibrium','correlated')
    if ( leqi(chars,'correlated') ) then
       TS_W_METHOD = TS_W_CORRELATED
    else if ( leqi(chars,'uncorrelated') ) then
       TS_W_METHOD = TS_W_UNCORRELATED
    else if ( leqi(chars,'k-uncorrelated') ) then
       TS_W_METHOD = TS_W_K_UNCORRELATED
    else
       call die('Could not determine flag TS.Weight.NonEquilibrium, please &
            &see manual.')
    end if

    ! Figure out the number of orbitals on the buffer atoms
    na_BufL     = fdf_get('TS.BufferAtomsLeft',na_BufL_def)
    no_BufL = 0
    do i = 1 , na_BufL
       no_BufL = no_BufL + lasto(i) - lasto(i-1)
    end do

    na_BufR = fdf_get('TS.BufferAtomsRight',na_BufR_def)
    no_BufR = 0
    do i = na_u - na_BufR + 1 , na_u
       no_BufR = no_BufR + lasto(i) - lasto(i-1)
    end do
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

    if ( TSmode ) then
       call read_contour_options( kT_in, IsVolt, VoltL, VoltR )
    end if



    call fdf_deprecated('TS.CalcGF','TS.ReUseGF')
    ReUseGF    = fdf_get('TS.ReUseGF',ReUseGF_def)
    UseVFix    = fdf_get('TS.UseVFix',UseVFix_def)
    ElecValenceBandBot = fdf_get('TS.Elec.Calc.BandBottom', &
         ElecValenceBandBot_def)

    ! To determine the same coordinate nature of the electrodes
    Elec_xa_EPS= fdf_get('TS.Elec.Coord.Eps',1e-4_dp,'Bohr')

    ! detect how many electrodes we have
    if ( fdf_nElec('TS',Elecs) > 0 ) then

       Elecs(:)%mu = 0._dp
       Elecs(1)%mu = VoltL
       Elecs(size(Elecs))%mu = VoltR

       do i = 1 , size(Elecs)
          ! Default things that could be of importance
          Elecs(i)%t_dir = ts_tdir
          Elecs(i)%UseBulk = UseBulk
          Elecs(i)%UpdateDMCR = UpdateDMCR
          if ( .not. fdf_Elec('TS',slabel,Elecs(i)) ) then
             call die('Could not find electrode: '//trim(name(Elecs(i))))
          end if
          ! set the placement in orbitals
          Elecs(i)%idx_no = lasto(Elecs(i)%idx_na-1)+1
       end do

       if ( sum(TotUsedAtoms(Elecs)) >= na_u ) then
          call die('Electrodes occupy the entire device')
       end if

       ! We need to sort the electrodes
       do i = 1 , size(Elecs) - 1
          idx = Elecs(i)%idx_na+TotUsedAtoms(Elecs(i))-1
          j = i + minloc(Elecs(i+1:)%idx_na,dim=1)
          if ( idx > Elecs(j)%idx_na ) then
             tmpElec = Elecs(j)
             Elecs(j) = Elecs(i)
             Elecs(i) = tmpElec
          end if
       end do

       ! we need to check that they indeed do not overlap
       do i = 1 , size(Elecs) - 1
          idx1 = Elecs(i)%idx_na
          idx2 = idx1 + TotUsedAtoms(Elecs(i)) - 1
          if ( idx1 <= na_BufL ) then
             print *,1,idx1,idx2,na_u
             call die('Buffer atoms overlap an electrode')
          else if ( na_u - na_BufR <= idx2 ) then
             print *,1,idx1,idx2,na_u
             call die('Buffer atoms overlap an electrode')
          end if
          do j = i + 1 , size(Elecs)
             ! if the index is smaller (then we have an error)
             if ( Elecs(j)%idx_na <= idx1 ) then
                call die('Sorting of electrodes went wrong, ensure no overlapping &
                     &electrodes')
             else if ( Elecs(j)%idx_na <= idx2 ) then
                call die('Overlapping electrodes is not physical, please correct.')
             end if
          end do
       end do
       ! check the last electrode
       i = size(Elecs)
       idx1 = Elecs(i)%idx_na
       idx2 = idx1 + TotUsedAtoms(Elecs(i)) - 1
       if ( idx1 <= na_BufL ) then
          print *,1,idx1,idx2,na_u
          call die('Buffer atoms overlap an electrode')
       else if ( na_u - na_BufR < idx2 ) then
          print *,1,idx1,idx2,na_u
          call die('Buffer atoms overlap an electrode')
       end if

       if ( size(Elecs) > 2 ) call die('currently does not work')

    else
       allocate(Elecs(2))
       ! defaults...
       Elecs(:)%t_dir      = ts_tdir
       Elecs(:)%UseBulk    = UseBulk
       Elecs(:)%UpdateDMCR = UpdateDMCR

       ! Setup the correct transmission directions
       Elecs(1)%inf_dir = INF_NEGATIVE
       Elecs(1)%Name    = 'Left'
       Elecs(1)%HSfile  = fdf_get('TS.HSFileLeft',HSFile_def)
       Elecs(1)%GFfile  = fdf_get('TS.GFFileLeft',trim(slabel)//'.TSGFL')
       Elecs(1)%GFtitle = fdf_get('TS.Elec.Left.GF.Title','Left Greens function')
       Elecs(1)%na_used = fdf_get('TS.NumUsedAtomsLeft',-1)
       call check_HSfile(Elecs(1))
       Elecs(1)%RepA1   = fdf_get('TS.ReplicateA1Left',1)
       Elecs(1)%RepA2   = fdf_get('TS.ReplicateA2Left',1)
       Elecs(1)%RepA3   = fdf_get('TS.ReplicateA3Left',1)
       if ( RepA1(Elecs(1)) < 1 .or. RepA2(Elecs(1)) < 1 .or. RepA3(Elecs(1)) < 1 ) &
            call die("Repetition in left electrode must be >= 1.")
       Elecs(1)%idx_na = na_BufL + 1
       Elecs(1)%idx_no = lasto(Elecs(1)%idx_na-1)+1
       Elecs(1)%mu = VoltL

       Elecs(2)%inf_dir = INF_POSITIVE
       Elecs(2)%Name    = 'Right'
       Elecs(2)%HSfile  = fdf_get('TS.HSFileRight',HSFile_def)
       Elecs(2)%GFfile  = fdf_get('TS.GFFileRight',trim(slabel)//'.TSGFR')
       Elecs(2)%GFtitle = fdf_get('TS.Elec.Right.GF.Title','Right Greens function')
       Elecs(2)%na_used = fdf_get('TS.NumUsedAtomsRight',-1)
       call check_HSfile(Elecs(2))
       Elecs(2)%RepA1   = fdf_get('TS.ReplicateA1Right',1)
       Elecs(2)%RepA2   = fdf_get('TS.ReplicateA2Right',1)
       Elecs(2)%RepA3   = fdf_get('TS.ReplicateA3Right',1)
       if ( RepA1(Elecs(2)) < 1 .or. RepA2(Elecs(2)) < 1 .or. RepA3(Elecs(2)) < 1 ) &
            call die("Repetition in left electrode must be >= 1.")
       Elecs(2)%idx_na  = na_u - na_BufR - TotUsedAtoms(Elecs(2)) + 1
       Elecs(2)%idx_no = lasto(Elecs(2)%idx_na-1)+1
       Elecs(2)%mu = VoltR

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
    
    ! Show the deprecated and obsolete labels
    call fdf_deprecated('TS.TriDiag','TS.SolutionMethod')
    call fdf_obsolete('TS.FixContactCharge')
    call fdf_obsolete('TS.KxyPoints')
    call fdf_obsolete('TS.NKVoltScale')

    ! Output Used Options in OUT file ....
    if (IOnode) then
       write(*,1) 'Save H and S matrices', saveTSHS
       write(*,1) 'Save S and quit (onlyS)', onlyS
    end if

    if (IONode .and. TSmode ) then
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
       write(*,7) 'Electronic temperature',kT/eV,'eV'
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
          write(*,6) 'Voltage', VoltFDF/eV,'Volts'
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
       write(*,1) 'Bulk Values in electrodes', UseBulk
       write(*,1) 'Update DM Contact Reg. only', UpdateDMCR

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
       write(*,1) 'Calc. band bottom in elec.', ElecValenceBandBot
       write(*,1) 'Re-use GF file if exists', ReUseGF
       write(*,10)'          >> Electrodes << '
       do i = 1 , size(Elecs)
          write(*,11)'>> '//trim(name(Elecs(i)))
          write(*,10)'  GF file', trim(GFFile(Elecs(i)))
          write(*,10)'  GF title', trim(GFtitle(Elecs(i)))
          write(*,10)'  Electrode TSHS file', trim(HSFile(Elecs(i)))
          write(*,5) '  # atoms used in electrode', UsedAtoms(Elecs(i))
          write(*,15)'  Electrode repetition A1/A2/A3', &
               RepA1(Elecs(i)),RepA2(Elecs(i)),RepA3(Elecs(i))
          if ( Elecs(i)%t_dir == 1 ) then
             chars = 'A1'
          else if ( Elecs(i)%t_dir == 2 ) then
             chars = 'A2'
          else if ( Elecs(i)%t_dir == 3 ) then
             chars = 'A3'
          end if
          write(*,5) '  Position in geometry', Elecs(i)%idx_na
          write(*,10) '  Transport direction for electrode', trim(chars)
          if ( Elecs(i)%inf_dir == INF_POSITIVE ) then
             write(*,10) '  Semi-infinite direction for electrode', 'positive'
          else
             write(*,10) '  Semi-infinite direction for electrode', 'negative'
          end if
          write(*,7) '  Chemical shift', Elecs(i)%mu/eV,'eV'
          write(*,1) '  Bulk values in electrode', Elecs(i)%UseBulk
          write(*,1) '  Update cross terms contact/electrode', .not. Elecs(i)%UpdateDMCR
       end do

       ! Print the contour information
       call ts_print_contour_options(cEq,cnEq, Eq_Eta, nEq_Eta,N_poles,IsVolt)
       
       write(*,11) repeat('*', 62)
       write(*,*)

       write(*,'(3a)') repeat('*',24),' Begin: TS CHECKS AND WARNINGS ',repeat('*',24)

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
       if ( IsVolt ) then
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

       ! print the warnings...
       call ts_print_contour_warnings(cEq,cnEq,kT, Eq_Eta, nEq_Eta, N_poles, IsVolt)

    end if

    if ( SaveTSHS .and. FixSpin ) then
       write(*,*) 'Fixed Spin not possible with Transiesta!'
       write(*,*) 'Electrodes with fixed spin is not possible with Transiesta !'
       call die('Stopping code')
    end if

    if ( IONode ) then
       write(*,'(3a,/)') repeat('*',24), &
            ' End: TS CHECKS AND WARNINGS ',repeat('*',26)
    end if

    if ( TSmode ) then

       ! Print out the contour blocks
       if ( IONode ) then
          write(*,'(/,a)') 'transiesta: contour input as perceived:'
       end if

       call ts_print_contour_block('TS.Contour.Eq',cEq, &
            msg_first='# First item is ALWAYS the circle part of the contour', &
            msg_last='# Last item is ALWAYS the tail part of the line contour')
       if ( IsVolt ) then
          if ( IONode ) write(*,*)
          call ts_print_contour_block('TS.Contour.nEq',cnEq, &
               msg_first='# First item is ALWAYS the lower tail part of the contour', &
               msg_last='# Last item is ALWAYS the upper tail part of the contour')
       end if
       if ( IONode ) write(*,'(/)')

    end if

    call fdf_obsolete('TS.GFTitle')

1   format('ts_options: ',a,t53,'=',4x,l1)
5   format('ts_options: ',a,t53,'=',i5,a)
6   format('ts_options: ',a,t53,'=',f10.4,tr1,a)
7   format('ts_options: ',a,t53,'=',f12.6,tr1,a)
8   format('ts_options: ',a,t53,'=',f10.4)
10  format('ts_options: ',a,t53,'=',4x,a)
11  format('ts_options: ',a)
15  format('ts_options: ',a,t53,'= ',i0,' X ',i0,' X ',i0)

contains 
  
  subroutine check_HSfile(El)
    type(Elec), intent(inout) :: el
    logical :: exist

    if ( TSmode ) then

       ! Check existance for left Electrode.TSHS
       inquire(file=TRIM(HSfile(El)),exist=exist)
       if ( .not. exist ) then
          call die("Electrode file does not exist. &
               &Please create electrode '"//trim(HSFile(El))//"' first.")
       end if
       ! Read in the number of atoms in the HSfile
       call ts_read_TSHS_opt(HSFile(El),no_u=El%no_u,na_u=El%na_u, &
            Bcast=.true.)
       allocate(El%xa(3,El%na_u),El%lasto(0:El%na_u))
       call ts_read_TSHS_opt(HSFile(El),xa=El%xa,lasto=El%lasto, &
            ucell=El%ucell,Ef=El%Ef, &
            Bcast=.true.)

       if ( UsedAtoms(El) < 0 ) then
          El%na_used = El%na_u
       else if ( UsedAtoms(El) == 0 ) then
          if(IONode) &
               write(*,*) "You need at least one atom in the electrode."
          call die("None atoms requested for electrode calculation.")
       else if ( El%na_u < UsedAtoms(El) ) then
          if (IONode) then
             write(*,*) "# of requested atoms is larger than available."
             write(*,*) "Requested: ",UsedAtoms(El)
             write(*,*) "Available: ",El%na_u
          end if
          call die("Error on requested atoms.")
       end if

       ! We have determined the number of atoms in the 
       ! TSHS file
       ! Read in lasto to determine the number of orbitals 
       allocate(El%xa_used(3,El%na_used),El%lasto_used(0:El%na_used))
       El%lasto_used(0) = 0
       El%no_used = 0
       if ( El%inf_dir == INF_NEGATIVE ) then ! old 'left'
          ! We use the last atoms
          do i = El%na_u - UsedAtoms(El) + 1 , El%na_u
             El%lasto_used(i) = El%lasto_used(i-1) + lasto(i)-lasto(i-1)
             El%xa_used(:,i)  = El%xa(:,i)
          end do
       else if ( El%inf_dir == INF_POSITIVE ) then ! old 'right'
          ! We use the last atoms
          do i = 1 , El%na_used
             El%lasto_used(i) = El%lasto_used(i-1) + lasto(i)-lasto(i-1)
             El%xa_used(:,i)  = El%xa(:,i)
          end do
       else
          call die('Unknown direction for the semi-infinite lead')
       end if
       El%no_used = El%lasto_used(El%na_used)

       ! We deallocate xa and lasto as they are not needed
       deallocate(El%xa,El%lasto)

    end if

  end subroutine check_HSfile

end subroutine read_ts_options

end module m_ts_options
