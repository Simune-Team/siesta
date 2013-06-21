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
real(dp) :: CCEmin       ! EMin for the Complex Contour (LB)
real(dp) :: GFEta        ! Imaginary part of the Bias Contour  
real(dp) :: kT           ! Electronic Temperature
integer  :: nline        ! Number of points on the "line" segment of the Contour
integer  :: ncircle      ! Number of points on the circle part of the contour 
integer  :: npol         ! Number of poles included in the contour
integer  :: nvolt        ! Number of points for the Bias integartion part
integer  :: ntransport   ! Number of points for transport calculation
integer  :: na_BufL      ! Number of Left Buffer Atoms
integer  :: na_BufR      ! Number of Right Buffer Atoms
integer  :: no_BufL      ! Number of Left Buffer orbitals
integer  :: no_BufR      ! Number of Right Buffer orbitals
character(200) :: GFTitle ! Title to paste in electrode Green's function files
character(200) :: GFFileL ! Electrode Left GF File
character(200) :: GFFileR ! Electrode Right GF File
type(Elec) :: ElLeft, ElRight
logical :: ElecValenceBandBot ! Calculate Electrode valence band bottom when creating electrode GF
integer :: Cmethod        ! Method for the contour integration
logical :: ReUseGF        ! Calculate the electrodes GF
logical :: ImmediateTSmode=.false. ! will determine to immediately start the transiesta
                           ! SCF. This is useful when you already have a converged
                           ! siesta DM

! If the energy-contour is not perfectly divisable by the number of nodes then adjust
logical :: NEn_Node_Correct = .true.
integer :: opt_TriMat_method = 0 ! Optimization method for determining the best tri-diagonal matrix split
! 0  == We optimize for speed
! 1  == We optimize for memory

logical :: VoltageInC ! Determines whether the voltage-drop should be located in the constriction
                      ! I.e. if the electrode starts at 10 Ang and the central region ends at 20 Ang
                      ! then the voltage drop will only take place between 10.125 Ang and 19.875 Ang

real(dp) :: Elec_xa_EPS

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
real(dp), parameter :: CCEmin_def = -3.0_dp  ! in Ry
real(dp), parameter :: GFEta_def = 0.000001_dp  ! in Ry
real(dp), parameter :: kT_def = 0.0019_dp  ! in Ry
integer, parameter :: nline_def = 6
integer, parameter :: ncircle_def = 24
integer, parameter :: npol_def = 6
integer, parameter :: nvolt_def = 5
integer, parameter :: ntransport_def = 0
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
  
  subroutine read_ts_options(ucell, na_u, lasto)

! SIESTA Modules Used
    use files, only  : slabel
    use fdf, only : fdf_get, fdf_deprecated, fdf_obsolete
    use fdf, only : leqi
    use parallel, only: IOnode, Nodes, operator(.parcount.)
    use units, only: eV
    use m_ts_contour, only : CC_METHOD_SOMMERFELD
    use m_ts_contour, only : CC_METHOD_GAUSSFERMI
    use m_ts_global_vars, only : ts_istep, TSinit
    use m_ts_io, only : ts_read_TSHS_na
    use m_ts_io, only : ts_read_TSHS_lasto
    use m_ts_method
    use m_ts_weight
    use m_ts_charge

    implicit none
    
    real(dp),intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, lasto(0:na_u)
! Internal Variables
    character(len=40) :: chars, s_cmethod
    integer :: i

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
    savetshs    = fdf_get('TS.SaveHS',savetshs_def)
    onlyS       = fdf_get('TS.onlyS',onlyS_def)

    VoltFDF     = fdf_get('TS.Voltage',voltfdf_def,'Ry') 
    IsVolt = dabs(VoltFDF) > 0.001_dp/eV

    ! Set up the fermi shifts for the left and right electrodes
    VoltL =  0.5_dp*VoltFDF
    VoltR = -0.5_dp*VoltFDF

    ! currently this does not work
    !ImmediateTSmode = fdf_get('TS.SCFImmediate',.false.)

    UseBulk     = fdf_get('TS.UseBulkInElectrodes',UseBulk_def)
    UpdateDMCR  = fdf_get('TS.UpdateDMCROnly',UpdateDMCR_def)
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
    chars = fdf_get('TS.Weight.NonEquilibrium','k-uncorrelated')
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
    TS_RHOCORR_FACTOR = fdf_get('TS.ChargeCorrectionFactor',.75_dp)
    if ( TS_RHOCORR_FACTOR < 0.0_dp .or. &
         1.0_dp < TS_RHOCORR_FACTOR) then
       call die("Charge correction factor must be in the range [0;1]")
    endif

    CCEMin = fdf_get('TS.ComplexContourEmin',CCEMin_def,'Ry')
    GFEta = fdf_get('TS.biasContour.Eta',GFEta_def,'Ry')
    if ( GFEta <= 0.d0) call die('ERROR: GFeta <= 0.0, we do not allow for &
         &using the advanced Greens function, please correct.')
    kT = fdf_get('ElectronicTemperature',kT_def,'Ry')

    s_cmethod = fdf_get('TS.biasContour.method',smethod_def)
    if ( leqi(s_cmethod,'sommerfeld') ) then
       Cmethod = CC_METHOD_SOMMERFELD
    else if ( leqi(s_cmethod,'gaussfermi') ) then
       Cmethod = CC_METHOD_GAUSSFERMI
    else
       Cmethod = 0 ! For producing error message later on
    end if
    npol       = fdf_get('TS.ComplexContour.NPoles',npol_def)
    ncircle    = fdf_get('TS.ComplexContour.NCircle',ncircle_def)
    nline      = fdf_get('TS.ComplexContour.NLine',nline_def)
    nvolt      = fdf_get('TS.biasContour.NumPoints',nvolt_def)
    NEn_Node_Correct = fdf_get('TS.Contour.NoEmptyCycles',.true.)
    if ( NEn_Node_Correct ) then
       i = npol+ncircle+nline
       ! We immediately correct the number of energy-points for the contour
       if ( mod(i,Nodes) /= 0 ) then
          ncircle = ncircle + Nodes - mod(i,Nodes)
       end if
       if ( mod(nvolt,Nodes) /= 0 ) then
          nvolt = nvolt + Nodes - mod(nvolt,Nodes)
       end if
    end if

    !Ntransport = fdf_get('TS.Contour.NTransport',Ntransport_def)
    GFTitle    = fdf_get('TS.GFTitle',GFTitle_def)
    chars = trim(slabel)//'.TSGFL'
    GFFileL    = fdf_get('TS.GFFileLeft',trim(chars))
    chars = trim(slabel)//'.TSGFR'
    GFFileR    = fdf_get('TS.GFFileRight',trim(chars))
    ReUseGF    = fdf_get('TS.ReUseGF',ReUseGF_def)
    UseVFix    = fdf_get('TS.UseVFix',UseVFix_def)
    ElecValenceBandBot = fdf_get('TS.CalcElectrodeValenceBandBottom', &
         ElecValenceBandBot_def)

    ! To determine the same coordinate nature of the electrodes
    Elec_xa_EPS= fdf_get('TS.Electrode.Coord.Eps',1e-4_dp,'Bohr')

    ElLeft%HSFile    = fdf_get('TS.HSFileLeft',HSFile_def)
    ElLeft%UsedAtoms = fdf_get('TS.NumUsedAtomsLeft',NUsedAtoms_def)
    call check_HSfile('Left',ElLeft)
    ElLeft%RepA1     = fdf_get('TS.ReplicateA1Left',NRepA_def)
    ElLeft%RepA2     = fdf_get('TS.ReplicateA2Left',NRepA_def)
    if ( RepA1(ElLeft) < 1 .or. RepA2(ElLeft) < 1 ) &
         call die("Repetition in left electrode must be >= 1.")
    
    ElRight%HSFile    = fdf_get('TS.HSFileRight',HSFile_def)
    ElRight%UsedAtoms = fdf_get('TS.NumUsedAtomsRight',NUsedAtoms_def)
    call check_HSfile('Right',ElRight)
    ElRight%RepA1     = fdf_get('TS.ReplicateA1Right',NRepA_def)
    ElRight%RepA2     = fdf_get('TS.ReplicateA2Right',NRepA_def)
    if ( RepA1(ElRight) < 1 .or. RepA2(ElRight) < 1 ) &
         call die("Repetition in right electrode must be >= 1.")


    ! Read in information about the voltage placement.
    chars = fdf_get('TS.VoltagePlacement','central')
    VoltageInC = .false.
    if ( leqi(trim(chars),'cell') ) then
       VoltageInC = .false.
    else if ( leqi(trim(chars),'central') .or. &
         leqi(trim(chars),'scat') ) then
       VoltageInC = .true.
    end if
    

    ! Here we check whether the user could perform the same
    ! calculation with the same GF-file
    ! We check that the user does not request the same GF files
    ! for runs with Bias. Furthermore, if na_u in Elec /= {NUsedAtomsL,NUsedAtomsR}
    ! then this is also not allowed.
    ! For non bias and na_u_elec == NUsedAtomsL == NUsedAtomsR
    ! then this is perfectly acceptable!
    if ( TSmode .and. trim(GFFileL) == trim(GFFileR) ) then ! Has to be case-sensitive !
       ! Read in the total number (if NumUsedAtoms is not the full...)
       call ts_read_TSHS_na(HSFile(ElLeft),i)

       ! They are the same
       if ( IsVolt ) call die("The same Green's function file &
            &can not be used in a bias calculation.")
       if ( trim(HSFile(ElLeft)) /= trim(HSFile(ElRight)) ) &
            call die("The same Green's function file &
            &can not be used if you request different &
            &electrode files.")
       if ( UsedAtoms(ElLeft) /= UsedAtoms(ElRight) .or. &
            UsedAtoms(ElLeft) /= i ) &
            call die("The same Green's function file &
            &can not be used if you do not request all &
            &atoms in the electrode!")
    else if ( (.not. IsVolt) .and. & ! for non-bias
         trim(HSFile(ElLeft)) == trim(HSFile(ElRight)) .and. & ! for same TSHS files
         UsedAtoms(ElLeft) == UsedAtoms(ElRight) .and. & ! for same number of atoms used
         i == UsedAtoms(ElLeft) ) then ! for using ALL atoms in the electrode
   ! For now this notification has been disabled, but in reality 
   ! could be enforced...
   !if ( IONode ) then
   !   write(*,*) 'NOTICE: In non-bias calculations, you can with'
   !   write(*,*) '        benefit use the same GF-files for both'
   !   write(*,*) '        the left and right electrode.'
   !   write(*,*) '        This *only* requires that you use ALL'
   !   write(*,*) '        atoms in the electrode and the left/right'
   !   write(*,*) '        TSHS files are the same.'
   !end if
    end if

    ! Show the deprecated and obsolete labels
    call fdf_deprecated('TS.CalcGF','TS.ReUseGF')
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
       write(*,1) 'Start TS-SCF cycle immediately', ImmediateTSmode
       if ( IsVolt ) then
          write(*,6) 'Voltage', VoltFDF/eV,' Volts'
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
       write(*,1) 'Bulk Values in Electrodes', UseBulk
       write(*,1) 'Update DM Contact Reg. only', UpdateDMCR

       write(*,5) 'Left buffer atoms', na_BufL
       write(*,5) 'Right buffer atoms', na_BufR

       write(*,5) 'N. Pts. Circle', ncircle
       write(*,5) 'N. Pts. Line', nline
       write(*,5) 'N. Poles in Contour', npol
       write(*,5) 'N. Pts. Bias Contour', nvolt
    ! write(*,5) 'N. Pts. Transport', Ntransport
       write(*,6) 'Contour E Min.', CCEmin / eV,' eV'

       write(*,7) 'GFEta', GFEta / eV,' eV'
       write(*,6) 'Electronic Temperature', kT / eV , ' eV'
       write(*,10)'Bias Contour Method', trim(s_cmethod)
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
       write(*,10)'GF title', trim(GFTitle)
       write(*,10)'Left GF File', trim(GFFileL)
       write(*,10)'Right GF File', trim(GFFileR)
       write(*,1) 'Re-use GF file if exists', ReUseGF
       write(*,10)'Left electrode TSHS file', trim(HSFile(ElLeft))
       write(*,5) '# atoms used in left elec. ', UsedAtoms(ElLeft)
       write(*,15) 'Left elec. repetition A1/A2', RepA1(ElLeft),RepA2(ElLeft)
       write(*,10)'Right electrode TSHS file', trim(HSFile(ElRight))
       write(*,5) '# atoms used in right elec. ', UsedAtoms(ElRight)
       write(*,15) 'Right elec. repetition A1/A2', RepA1(ElRight),RepA2(ElRight)


       write(*,11) repeat('*', 62)
       write(*,*)

       write(*,'(3a)') repeat('*',24),' Begin: TS CHECKS AND WARNINGS ',repeat('*',24)

       ! Check that the unitcell does not extend into the transport direction
       do i = 1 , 2
          if ( abs(ucell(3,i)) > 1e-7 .or. abs(ucell(i,3)) > 1e-7 ) then
             write(*,*) &
                  "ERROR: Unitcell has the electrode extend into the &
                  &transport direction."
             write(*,*) &
                  "Please change the geometry."
             call die("Electrodes extend into the transport direction. &
                  &Please change the geometry.")
          end if
       end do

! Print out message if the number of contour points are not 
! divisable by the number of Nodes
! We have two cases
! IsVolt:
!   - The equilibrium parts are the same computational cost
!   - The voltage contour point is more "heavy" in computation
!   * Solution make both Left and Right equi contours divisible by Nodes
!   * Make Nvolt divisible by Nodes
       if ( IsVolt ) then
          if ( mod(2*(Npol+NLine+Ncircle),Nodes) /= 0 ) then
             write(*,*) "NOTICE: Equilibrium energy contour points are not"
             write(*,*) "        divisable by the number of nodes."
             write(*,*) "        Better scalability is achived by changing:"
             write(*,*) "          - TS.ComplexContour.NPoles"
             write(*,*) "          - TS.ComplexContour.NLine"
             write(*,*) "          - TS.ComplexContour.NCircle"

             ! Calculate optimal number of energy points
             i = 2*(Npol+Nline+Ncircle)
             write(*,'(t10,a,i4)') "Used equilibrium # of energy points   : ",i
             i = Nodes .PARCOUNT. i
             write(*,'(t10,a,i4,tr1,a4,i3,/)') &
                  "Optimal equilibrium # of energy points: ",i, &
                  "+- i*",Nodes
          end if
          if ( mod(NVolt,Nodes) /= 0 ) then
             write(*,*) "NOTICE: Non-equilibrium energy contour points are not"
             write(*,*) "        divisable by the number of nodes."
             write(*,*) "        Better scalability is achieved by changing:"
             write(*,*) "          - TS.ComplexContour.NVolt"
          
             ! Calculate optimal number of energy points
             i = NVolt
             write(*,'(t10,a,i4)') "Used non-equilibrium # of energy points   : ",i
             i = Nodes .PARCOUNT. i
             write(*,'(t10,a,i4,tr1,a4,i3,/)') &
                  "Optimal non-equilibrium # of energy points: ",i, &
                  "+- i*",Nodes
          end if
          if ( mod(2*(Npol+Nline+Ncircle)+NVolt,Nodes) /= 0 ) then
             write(*,*) "NOTICE: Total energy contour points are not"
             write(*,*) "        divisable by the number of nodes."

             ! Calculate optimal number of energy points
             i = 2*(Npol+Nline+Ncircle)+NVolt
             write(*,'(t10,a,i4)') "Used # of energy points   : ",i
             i = Nodes .PARCOUNT. i
             write(*,'(t10,a,i4,tr1,a4,i3,/)') &
                  "Optimal # of energy points: ",i,"+- i*",Nodes
          end if
       else
! .not. IsVolt:
!   - The equilibrium parts are the same computational cost
!   * Solution make the equi contours divisible by Nodes
          if ( mod(Npol+NLine+Ncircle,Nodes) /= 0 ) then
             write(*,*) "NOTICE: Total number of energy points is &
                  &not divisable by the number of nodes."
             write(*,*) "        There are no computational costs &
                  &associated with increasing this."
! Calculate optimal number of energy points
             i = Npol+Nline+Ncircle
             write(*,'(t10,a,i4)') "Used # of energy points   : ",i
             i = Nodes .PARCOUNT. i
             write(*,'(t10,a,i4)') "Optimal # of energy points: ",i
          end if
       end if

    end if

    ! Integration Method
    if( Cmethod == 0 ) then
       if ( IONode ) then
          write(*,*) 'WARNING: TS.biasContour.method not recognized.'
          write(*,*) '         Reverting to gaussfermi instead'
       end if
       Cmethod = CC_METHOD_GAUSSFERMI
    end if

    if (fixspin ) then
       write(*,*) 'Fixed Spin not possible in TS Calculations !'
       call die('Stopping code')
    end if

    if ( IONode ) then
       write(*,'(3a,/)') repeat('*',24), &
            ' End: TS CHECKS AND WARNINGS ',repeat('*',26)
    end if

1   format('ts_options: ',a,t53,'=',4x,l1)
5   format('ts_options: ',a,t53,'=',i5,a)
6   format('ts_options: ',a,t53,'=',f10.4,a)
7   format('ts_options: ',a,t53,'=',f12.6,a)
8   format('ts_options: ',a,t53,'=',f10.4)
10  format('ts_options: ',a,t53,'=',4x,a)
11  format('ts_options: ',a)
15  format('ts_options: ',a,t53,'= ',i0,' X ',i0)

contains 
  
  subroutine check_HSfile(LR,el)
    character(len=*), intent(in) :: LR
    type(Elec), intent(inout) :: el
    integer :: tmp_NUsedAtoms
    integer, allocatable, dimension(:) :: lasto
    logical :: exist
    if ( TSmode ) then
! Check existance for left Electrode.TSHS
       inquire(file=TRIM(HSfile(el)),exist=exist)
       if ( .not. exist ) then
          call die(trim(LR)//" electrode file does not exist. &
               &Please create electrode '"//trim(HSFile(el))//"' first.")
       end if
       ! Read in the number of atoms in the HSfile
       call ts_read_TSHS_na(HSFile(el),tmp_NUsedAtoms)

       if ( UsedAtoms(el) < 0 ) then
          el%UsedAtoms = tmp_NUsedAtoms
       else if ( UsedAtoms(el) == 0 ) then
          if(IONode) &
               write(*,*) "You need at least one atom in the electrode."
          call die("None atoms requested for electrode calculation.")
       else if ( tmp_NUsedAtoms < UsedAtoms(el) ) then
          if (IONode) then
             write(*,*) "# of requested atoms is larger than available."
             write(*,*) "Requested: ",UsedAtoms(el)
             write(*,*) "Available: ",tmp_NUsedAtoms
          end if
          call die("Error on requested atoms.")
       end if

       ! We have determined the number of atoms in the 
       ! TSHS file
       ! Read in lasto to determine the number of orbitals 
       ! used in the electrode
       allocate(lasto(0:tmp_NUsedAtoms))
       call ts_read_TSHS_lasto(HSFile(el),tmp_NUsedAtoms,lasto)
       el%UsedOrbs = 0
       if ( LR == 'Left' ) then
          ! We use the first atoms
          do i = 1 , UsedAtoms(el)
             el%UsedOrbs = el%UsedOrbs + lasto(i)-lasto(i-1)
          end do
       else
          ! We use the last atoms
          do i = tmp_NUsedAtoms - UsedAtoms(el) + 1 , tmp_NUsedAtoms
             el%UsedOrbs = el%UsedOrbs + lasto(i)-lasto(i-1)
          end do
       end if
    end if

  end subroutine check_HSfile
  
end subroutine read_ts_options

end module m_ts_options
