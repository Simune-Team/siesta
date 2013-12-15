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
  
logical  :: SaveTSHS     ! Saves the Hamiltonian and Overlap matrices if the 
                         ! the option TS.SaveHS is specified in the input file
logical  :: onlyS        ! Option to only save overlap matrix
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
logical  :: IsVolt       ! Logical for dabs(VoltFDF) > 0.0001d*eV
real(dp) :: Volt         ! Bias applied, Internally Volt=voltfdf/eV (eV). 
real(dp) :: VoltL        ! Bias on the left electrode   (  .5 * Volt )
real(dp) :: VoltR        ! Bias on the right electrode  ( -.5 * Volt )
integer  :: na_BufL      ! Number of Left Buffer Atoms
integer  :: na_BufR      ! Number of Right Buffer Atoms
integer  :: no_BufL      ! Number of Left Buffer orbitals
integer  :: no_BufR      ! Number of Right Buffer orbitals
type(Elec), allocatable, target :: Elecs(:)
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

real(dp) :: Elecs_xa_EPS

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
  
  subroutine read_ts_options( wmix, kT, ucell, na_u, xa, lasto)

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
    
    real(dp), intent(in) :: wmix, kT
    real(dp),intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
! Internal Variables
    real(dp) :: tmp
    character(len=70) :: c, chars
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
    saveTSHS = fdf_get('TS.SaveHS',.true.)
    onlyS    = fdf_get('TS.onlyS',.false.)

    Volt     = fdf_get('TS.Voltage',0._dp,'Ry') 
    ! Voltage situation is above 0.1 meV (probably too low...)
    IsVolt   = dabs(Volt) > 0.0001_dp*eV
    if ( .not. IsVolt ) Volt = 0._dp

    ! Set up the fermi shifts for the left and right electrodes
    ! in case of several electrodes the first and last electrode
    ! correspond to left/right respectively. (the remaining electrodes have mu=0)
    VoltL =  0.5_dp * Volt
    VoltR = -0.5_dp * Volt

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
    na_BufL     = fdf_get('TS.BufferAtomsLeft',0)
    no_BufL = 0
    do i = 1 , na_BufL
       no_BufL = no_BufL + lasto(i) - lasto(i-1)
    end do

    na_BufR = fdf_get('TS.BufferAtomsRight',0)
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

    call fdf_deprecated('TS.CalcGF','TS.ReUseGF')
    ReUseGF    = fdf_get('TS.ReUseGF',.false.)

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
       call fdf_obsolete('TS.ReplicateA3'//trim(chars))
    end do

    ! notice that this does not have the same meaning... 
    call fdf_deprecated('TS.UpdateDMCROnly','TS.Elecs.DM.CrossTerms')
    call fdf_deprecated('TS.UseBulk','TS.Elecs.Bulk')

    ! To determine the same coordinate nature of the electrodes
    Elecs_xa_EPS= fdf_get('TS.Elecs.Coord.Eps',1e-4_dp,'Bohr')

    ! detect how many electrodes we have
    if ( fdf_nElec('TS',Elecs) < 2 ) then
       ! The user has done something wrong...
       ! Write out a working FDF (similar to the first thing)
       write(*,*) '%block TS.Elecs'
       write(*,*) '  Left'
       write(*,*) '  Right'
       write(*,*) '%endblock TS.Elecs'
       do i = 1 , 2
          if ( i == 1 ) then
             chars = 'Left'
             c = 'negative'
             j = na_BufL + 1
             tmp = VoltL
          else
             chars = 'Right'
             c = 'positive'
             j = na_u - na_BufR
             tmp = VoltR
          end if
          
          write(*,*)
          write(*,*) '%block TS.Elec.'//trim(chars)
          write(*,*) '  TSHS <TSHS-file for '//trim(chars)//' electrode> ('//trim(chars)//'.TSHS)'
          write(*,*) '  semi-inf-direction <direction of infinity for '//trim(chars)//' electrode> ('//trim(c)//')'
          write(*,*) '  chemical-shift <chemical potential in '//trim(chars)//' electrode> (',tmp/eV,' eV)'
          write(*,*) '  electrode-position <position in FDF structure of '//trim(chars)//' electrode> (',j,')'
          write(*,*) '  #   Other options (see the manual if in doubt):'
          write(*,*) '  # bulk <Electrode is bulk>'
          write(*,*) '  # cross-terms <whether the DM for cross-terms gets updated>'
          write(*,*) '  # GF-title <title for GF file>'
          write(*,*) '  # GF <file name for GF>'
          write(*,*) '  # used-atoms <number of atoms used in TSHS file>'
          write(*,*) '  # replicate-a1 <repetition in a1-direction>'
          write(*,*) '  # replicate-a2 <repetition in a2-direction>'
          write(*,*) '  # replicate-a3 <repetition in a3-direction>'
          write(*,*) '  # replicate <rep in a1> <rep in a2> <rep in a3>'
          write(*,*) '  contour.eq'
          write(*,*) '    begin'
          write(*,*) '      C-'//trim(chars)
          write(*,*) '      L-'//trim(chars)
          write(*,*) '      T-'//trim(chars)
          write(*,*) '    end'
          !write(*,*) '  # transport (please see the manual)'
          write(*,*) '%endblock TS.Elec.'//trim(chars)
       end do
       write(*,*)
       
       write(*,*) ' *** CREATE DEFAULT CONTOUR BLOCK **** '

       call die('Please see the output for how to construct an example electrode configuration')
    end if

    ! Setup default parameters for the electrodes
    ! first electrode is the "left"
    ! last electrode is the "right"
    ! the remaining electrodes have their chemical potential at 0
    ! Currently the transport direction for all electrodes is the default
    ! We should probably warn if +2 electrodes are used and t_dir is the
    ! same for all electrodes... Then the user needs to know what (s)he is doing...
    Elecs(:)%mu = 0._dp
    Elecs(1)%mu = VoltL
    Elecs(size(Elecs))%mu = VoltR
    Elecs(:)%t_dir = ts_tdir
    Elecs(:)%Bulk = fdf_get('TS.Elecs.Bulk',.true.)
    if ( .not. Elecs(1)%Bulk ) then
       Elecs(:)%DM_CrossTerms = .true.
    else
       Elecs(:)%DM_CrossTerms = fdf_get('TS.Elecs.DM.CrossTerms',.false.)
    end if
    ! We default to not calculate the band-bottom...
    Elecs(:)%BandBottom = fdf_get('TS.Elecs.Calc.BandBottom', .false.)

    do i = 1 , size(Elecs)
       ! Default things that could be of importance
       if ( .not. fdf_Elec('TS',slabel,Elecs(i)) ) then
          call die('Could not find electrode: '//trim(name(Elecs(i))))
       end if
       ! set the placement in orbitals
       if ( Elecs(i)%idx_na < 0 ) &
            Elecs(i)%idx_na = na_u + Elecs(i)%idx_na
       if ( Elecs(i)%idx_na < 1 .or. &
            na_u < Elecs(i)%idx_na ) &
            call die("Electrode position does not exist")
       Elecs(i)%idx_no = lasto(Elecs(i)%idx_na-1)+1

       ! Populate the intrinsic things
       call check_HSfile(Elecs(i))
       
    end do

    ! check that all have at least 2 contour points on the equilibrium contour
    if ( .not. all(Eq_segs(Elecs) > 1) ) then
       print *,Eq_segs(Elecs)
       call die('All electrodes does not have at least 2 equilibrium contours')
    end if
    
    if ( sum(TotUsedAtoms(Elecs)) >= na_u ) then
       call die('Electrodes occupy the entire device')
    end if
       
    ! Sort the electrodes according to the structure
    do i = 1 , size(Elecs) - 1
       idx = Elecs(i)%idx_na+TotUsedAtoms(Elecs(i))-1
       j = i + minloc(Elecs(i+1:)%idx_na,dim=1)
       if ( idx > Elecs(j)%idx_na ) then
          tmpElec  = Elecs(j)
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
             print *,Elecs(j)%idx_na, idx1
             call die('Sorting of electrodes went wrong, ensure no overlapping &
                  &electrodes')
          else if ( Elecs(j)%idx_na <= idx2 ) then
             print *,Elecs(j)%idx_na, idx2
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

    ! CHECK THIS (we could allow it by only checking the difference...)
    write(*,*) 'TODO MU'
    if ( maxval(Elecs(:)%mu) > max(VoltL,VoltR) ) then
       call die('MU')
    else if ( minval(Elecs(:)%mu) < min(VoltL,VoltR) ) then
       call die('MU')
    end if

    ! WILL WORK EVENTUALLY
    write(*,*) 'TODO SEVERAL ELECTRODES'
    !if ( size(Elecs) > 2 ) call die('currently does not work')


    ! read in contour options
    if ( TSmode ) then
       call read_contour_options( Elecs, kT, IsVolt, Volt )
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
          write(*,7)  '  Chemical shift', Elecs(i)%mu/eV,'eV'
          write(*,1)  '  Bulk values in electrode', Elecs(i)%Bulk
          write(*,1)  '  Update cross terms contact/electrode', Elecs(i)%DM_CrossTerms
          write(*,1)  '  Calc. valence band-bottom eigenvalue', Elecs(i)%BandBottom
       end do

       ! Print the contour information
       call print_contour_options( 'TS' , IsVolt )
       
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
       call print_contour_block( 'TS' , IsVolt )

       ! write out the contour
       call io_contour(IsVolt, Elecs, slabel)

       ! Print out the electrode coordinates
       do i = 1 , size(Elecs)
          call print_elec(Elecs(i),na_u,xa)
       end do

    end if

1   format('ts_options: ',a,t53,'=',4x,l1)
5   format('ts_options: ',a,t53,'=',i5,a)
20  format('ts_options: ',a,t53,'= ',i0,' -- ',i0)
6   format('ts_options: ',a,t53,'=',f10.4,tr1,a)
7   format('ts_options: ',a,t53,'=',f12.6,tr1,a)
8   format('ts_options: ',a,t53,'=',f10.4)
10  format('ts_options: ',a,t53,'=',4x,a)
11  format('ts_options: ',a)
15  format('ts_options: ',a,t53,'= ',i0,' x ',i0,' x ',i0)

contains 
  
  subroutine check_HSfile(El)
    type(Elec), intent(inout) :: el
    logical :: exist
    integer :: i

    if ( TSmode ) then

       ! Check existance for left Electrode.TSHS
       inquire(file=TRIM(HSfile(El)),exist=exist)
       if ( .not. exist ) then
          call die("Electrode file does not exist. &
               &Please create electrode '"//trim(HSFile(El))//"' first.")
       end if

       ! Read in the number of atoms in the HSfile
       call ts_read_TSHS_opt(HSFile(El),no_u=El%no_u,na_u=El%na_u, &
            nspin=El%nspin, Ef=El%Ef, ucell=El%ucell, &
            Bcast=.true.)
       allocate(El%xa(3,El%na_u),El%lasto(0:El%na_u))
       call ts_read_TSHS_opt(HSFile(El),xa=El%xa,lasto=El%lasto, &
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
