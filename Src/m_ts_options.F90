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
logical  :: TriDiag      ! true if tridiagonalization
logical  :: UpdateDMCR   ! Update DM values of ONLY Central Region
integer  :: ChargeCorr   ! Integer holding the method of charge correction
                         !  0 => Will not do charge correction
                         !  1 => Excess/missing charge is corrected in the buffer layers
real(dp) :: ChargeCorr_factor ! A factor for the correction (should be in the range 0 <= 1)
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
integer  :: NBufAtL      ! Number of Left Buffer Atoms
integer  :: NBufAtR      ! Number of Right Buffer Atoms
integer  :: NRepA1L      ! Number of Left Repetitions in A1 direction
integer  :: NRepA2L      ! Number of Left Repetitions in A2 direction
integer  :: NRepA1R      ! Number of Right Repetitions in A1 direction
integer  :: NRepA2R      ! Number of Right Repetitions in A2 direction
integer  :: NUsedAtomsL  ! Number of atoms used from the Left electrode
integer  :: NUsedAtomsR  ! Number of atoms used from the Right electrode
integer  :: NUsedOrbsL   ! Number of orbitals used from the Left electrode
integer  :: NUsedOrbsR   ! Number of orbitals used from the Right electrode
character(200) :: GFTitle ! Title to paste in electrode Green's function files
character(200) :: GFFileL ! Electrode Left GF File
character(200) :: GFFileR ! Electrode Right GF File
character(200) :: HSFileL ! Electrode Left TSHS File
character(200) :: HSFileR ! Electrode Right TSHS File
logical       :: ElecValenceBandBot ! Calculate Electrode valence band bottom when creating electrode GF
integer :: Cmethod        ! Method for the contour integration
logical :: ReUseGF        ! Calculate the electrodes GF
logical :: ImmediateTSmode=.false. ! will determine to immediately start the transiesta
                           ! SCF. This is useful when you already have a converged
                           ! siesta DM

logical :: VoltageInC ! Determines whether the voltage-drop should be located in the constriction
                      ! I.e. if the electrode starts at 10 Ang and the central region ends at 20 Ang
                      ! then the voltage drop will only take place between 10.125 Ang and 19.875 Ang

!==========================================================================*
!==========================================================================*
!  Default Values for arguments read from input file                       *
!--------------------------------------------------------------------------*

logical, parameter :: savetshs_def = .true.
logical, parameter :: onlyS_def = .false.
logical, parameter :: tsdme_def = .true.
logical, parameter :: UseBulk_def = .true.
logical, parameter :: TriDiag_def = .true.
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
integer, parameter :: NBufAtL_def = 0
integer, parameter :: NBufAtR_def = 0
integer, parameter :: NRepA_def = 1
integer, parameter :: NUsedAtoms_def = -1
character(20), parameter :: smethod_def = 'gaussfermi'
character(33), parameter :: GFTitle_def = 'Generated GF file'
character(33), parameter :: HSFile_def = 'NOT REQUESTED'
character(4),  parameter :: ChargeCorr_def = 'none'
real(dp),  parameter :: ChargeCorr_factor_def = 0.75_dp
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
!
! **************************** OUTPUT *********************************

subroutine read_ts_options(ucell)

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
#ifdef MPI
use mpi_siesta, only : MPI_Integer, MPI_Comm_World
#endif
implicit none

real(dp),intent(in) :: ucell(3,3)
! Internal Variables
character(len=40) :: chars, s_cmethod
integer :: i
#ifdef MPI
integer :: MPIerror
#endif

if (isolve.eq.SOLVE_TRANSI) then
   TSmode = .true.
   ! If in TSmode default to initalization
   ! In case of 'DM.UseSaveDM TRUE' TSinit will be set accordingly
   TSinit = .true.
endif

if (IOnode) then
 write(*,'(/,2a)') 'ts_read_options: ', repeat('*', 62)
end if

!Set ts_istep default
ts_istep = 0

! Reading the Transiesta solution method
chars = fdf_get('TS.SolutionMethod','sparse')
if ( leqi(chars,'original') ) then
   ts_method = TS_ORIGINAL
else if ( leqi(chars,'sparse') ) then
   ts_method = TS_SPARSITY
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
TriDiag     = fdf_get('TS.TriDiag',TriDiag_def)
UpdateDMCR  = fdf_get('TS.UpdateDMCROnly',UpdateDMCR_def)
NBufAtL     = fdf_get('TS.BufferAtomsLeft',NBufAtL_def)
NBufAtR     = fdf_get('TS.BufferAtomsRight',NBufAtR_def)
if ( NBufAtL < 0 .or. NBufAtR < 0 ) then
   call die("Buffer atoms must be 0 or a positive integer.")
end if
chars       = fdf_get('TS.ChargeCorrection',ChargeCorr_def)
ChargeCorr = 0
if ( leqi(chars,'none') ) then
   ChargeCorr = 0
else if ( leqi(chars,'b') .or. leqi(chars,'buffer') ) then
   ChargeCorr = 1
end if
ChargeCorr_factor = fdf_get('TS.ChargeCorrectionFactor',ChargeCorr_factor_def)
if ( ChargeCorr_factor < 0.0_dp .or. &
     1.0_dp < ChargeCorr_factor) then
   call die("Charge correction factor must be in the range [0;1]")
endif
CCEMin     = fdf_get('TS.ComplexContourEmin',CCEMin_def,'Ry')
GFEta      = fdf_get('TS.biasContour.Eta',GFEta_def,'Ry')
if ( GFEta <= 0.d0) call die('ERROR: GFeta <= 0.0 ')
kT         = fdf_get('ElectronicTemperature',kT_def,'Ry')
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

HSFileL     = fdf_get('TS.HSFileLeft',HSFile_def)
NUsedAtomsL = fdf_get('TS.NumUsedAtomsLeft',NUsedAtoms_def)
call check_HSfile('Left',HSFileL,NUsedAtomsL,NUsedOrbsL)
NRepA1L     = fdf_get('TS.ReplicateA1Left',NRepA_def)
NRepA2L     = fdf_get('TS.ReplicateA2Left',NRepA_def)
if ( NRepA1L < 1 .or. NRepA2L < 1 ) &
     call die("Repetition in left electrode must be >= 1.")

HSFileR     = fdf_get('TS.HSFileRight',HSFile_def)
NUsedAtomsR = fdf_get('TS.NumUsedAtomsRight',NUsedAtoms_def)
call check_HSfile('Right',HSFileR,NUsedAtomsR,NUsedOrbsR)
NRepA1R     = fdf_get('TS.ReplicateA1Right',NRepA_def)
NRepA2R     = fdf_get('TS.ReplicateA2Right',NRepA_def)
if ( NRepA1R < 1 .or. NRepA2R < 1 ) &
     call die("Repetition in right electrode must be >= 1.")

chars = fdf_get('TS.VoltagePlacement','central')
VoltageInC = .false.
if ( leqi(trim(chars),'cell') ) then
   VoltageInC = .false.
else if ( leqi(trim(chars),'central') .or. &
     leqi(trim(chars),'scat') ) then
   VoltageInC = .true.
end if


! Setup the correct handling of EQUILIBRIUM solution method:
! See above the global variable for its use.
GF_INV_EQUI_PART = UseBulk .and. UpdateDMCR

! Here we check whether the user could perform the same
! calculation with the same GF-file
! We check that the user does not request the same GF files
! for runs with Bias. Furthermore, if na_u in Elec /= {NUsedAtomsL,NUsedAtomsR}
! then this is also not allowed.
! For non bias and na_u_elec == NUsedAtomsL == NUsedAtomsR
! then this is perfectly acceptable!
if ( TSmode .and. trim(GFFileL) == trim(GFFileR) ) then ! Has to be case-sensitive !
   ! Read in the total number (if NumUsedAtoms is not the full...)
   call ts_read_TSHS_na(HSFileL,i)

   ! They are the same
   if ( IsVolt ) call die("The same Green's function file &
        &can not be used in a bias calculation.")
   if ( trim(HSFileL) /= trim(HSFileR) ) &
        call die("The same Green's function file &
        &can not be used if you request different &
        &electrode files.")
   if ( NUsedAtomsL /= NUsedAtomsR .or. &
        NUsedAtomsL /= i ) &
        call die("The same Green's function file &
        &can not be used if you do not request all &
        &atoms in the electrode!")
else if ( (.not. IsVolt) .and. & ! for non-bias
     trim(HSFileL) == trim(HSFileR) .and. & ! for same TSHS files
     NUsedAtomsL == NUsedAtomsR .and. & ! for same number of atoms used
     i == NUsedAtomsL ) then ! for using ALL atoms in the electrode
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
call fdf_obsolete('TS.FixContactCharge')
call fdf_obsolete('TS.KxyPoints')
call fdf_obsolete('TS.NKVoltScale')


! Output Used Options in OUT file ....
if (IOnode) then
if ( ts_method == TS_ORIGINAL ) then
 write(*,10)'ts_read_options: Solution method              =', 'Original'
else if ( ts_method == TS_SPARSITY ) then
 write(*,10)'ts_read_options: Solution method              =', 'Sparsity pattern'
end if
 write(*,1) 'ts_read_options: Start SCF cycle immediately  =', ImmediateTSmode
 write(*,1) 'ts_read_options: Save H and S matrices        =', saveTSHS
 write(*,1) 'ts_read_options: Save S and quit (onlyS)      =', onlyS
end if
if (ionode .and. TSmode ) then
if ( IsVolt ) then
 write(*,6) 'ts_read_options: TranSIESTA Voltage           =', VoltFDF/eV,' Volts'
 if ( VoltageInC ) then
    write(*,'(a)')'ts_read_options: Voltage drop across central region'
 else
    write(*,'(a)')'ts_read_options: Voltage drop across entire cell'    
 end if
else
 write(*,'(a)')'ts_read_options: TranSIESTA no voltage applied'
end if
 write(*,1) 'ts_read_options: Bulk Values in Electrodes    =', UseBulk
 write(*,1) 'ts_read_options: TriDiag                      =', TriDiag 
 write(*,1) 'ts_read_options: Update DM Contact Reg. only  =', UpdateDMCR
 write(*,5) 'ts_read_options: N. Buffer At. Left           =', NBufAtL
 write(*,5) 'ts_read_options: N. Buffer At. Right          =', NBufAtR
 write(*,5) 'ts_read_options: N. Pts. Circle               =', ncircle
 write(*,5) 'ts_read_options: N. Pts. Line                 =', nline
 write(*,5) 'ts_read_options: N. Poles in Contour          =', npol
 write(*,5) 'ts_read_options: N. Pts. Bias Contour         =', nvolt
! write(*,5) 'ts_read_options: N. Pts. Transport            =', Ntransport
 write(*,6) 'ts_read_options: Contour E Min.               =', CCEmin,' Ry'
 write(*,7) 'ts_read_options: GFEta                        =', GFEta,' Ry'
 write(*,6) 'ts_read_options: Electronic Temperature       =', kT, ' Ry'
 write(*,10)'ts_read_options: Bias Contour Method          =', trim(s_cmethod)
if ( ChargeCorr == 0 ) then
 write(*,'(a)')'ts_read_options: Will not correct charge fluctuations'
else if ( ChargeCorr == 1 ) then ! Correct in buffer
  if ( 0 < NBufAtL .or. 0 < NBufAtR ) then
   write(*,10)'ts_read_options: Charge fluctuation correction=','buffer'
  else
     call die('Charge correction can not happen in buffer as no buffer &
          &atoms exist.')
  end if
  write(*,8)'ts_read_options: Charge correction factor     =',ChargeCorr_factor
end if
  write(*,1) 'ts_read_options: Calc. band bottom in elec.   =', ElecValenceBandBot
  write(*,10)'ts_read_options: GF title                     =', trim(GFTitle)
  write(*,10)'ts_read_options: Left GF File                 =', trim(GFFileL)
  write(*,10)'ts_read_options: Right GF File                =', trim(GFFileR)
  write(*,1) 'ts_read_options: Re-use GF file if exists     =', ReUseGF
  write(*,10)'ts_read_options: Left electrode TSHS file     =', trim(HSFileL)
  write(*,5) 'ts_read_options: # atoms used in left elec.   = ', NUsedAtomsL
  write(*,'(a,i3,'' X '',i3)') &
             'ts_read_options: Left elec. repetition A1/A2  = ', NRepA1L,NRepA2L
  write(*,10)'ts_read_options: Right electrode TSHS file    =', trim(HSFileR)
  write(*,5) 'ts_read_options: # atoms used in right elec.  = ', NUsedAtomsR
  write(*,'(a,i3,'' X '',i3)') &
             'ts_read_options: Right elec. repetition A1/A2 = ', NRepA1R,NRepA2R


 write(*,'(2a,/)') 'ts_read_options: ', repeat('*', 62)

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

! UseBulk and TriDiag
if ( ts_method == TS_ORIGINAL ) then
   if ( (.not. UseBulk) .and. TriDiag ) then
      if ( IONode ) then
         write(*,*) "WARNING: TriDiag only for UseBulkInElectrodes"
         write(*,*) "         Reverting to normal inversion scheme"
      end if
      TriDiag = .false. 
   end if
   if ( TriDiag .and. IsVolt .and. (.not. UpdateDMCR) ) then
      if ( IONode ) then
         write(*,*) "WARNING: TriDiag does not perform correctly in the"
         write(*,*) "         original solution method and with a bias."
         write(*,*) "         Reverting to normal inversion scheme"
      end if
      TriDiag = .false.
   end if
end if

! sparsity pattern
if ( ts_method == TS_SPARSITY ) then
   if ( TriDiag ) then
      ! Change to the correct method
      ts_method = TS_SPARSITY_TRI
   end if
   if ( ChargeCorr /= 0 ) then
      if ( IONode ) then
       write(*,*) "WARNING: Charge correction only for original solution method"
       write(*,*) "         No correction will be performed"
    end if
    ChargeCorr = 0
   end if
end if

! Integration Method
if( Cmethod == 0 ) then
   if ( IONode ) then
    write(*,*) 'WARNING: TS.biasContour.method not recognized.'
    write(*,*) '         Reverting to gaussfermi instead'
   end if
   Cmethod = CC_METHOD_GAUSSFERMI
endif

if (fixspin ) then
   write(*,*) 'Fixed Spin not possible in TS Calculations !'
   call die('Stopping code')
end if

if ( IONode ) then
   write(*,'(3a,/)') repeat('*',24), &
        ' End: TS CHECKS AND WARNINGS ',repeat('*',26)
end if

1   format(a,4x,l1)
5   format(a,i5,a)
6   format(a,f10.4,a)
7   format(a,f12.6,a)
8   format(a,f10.4)
10  format(a,4x,a)

contains 
  
  subroutine check_HSfile(LR,HSFile,NUsedAtoms,NUsedOrbs)
    character(len=*), intent(in) :: LR
    character(len=*), intent(in) :: HSFile
    integer, intent(inout) :: NUsedAtoms, NUsedOrbs
    integer :: tmp_NUsedAtoms
    integer, allocatable, dimension(:) :: lasto
    logical :: exist
    if ( TSmode ) then
! Check existance for left Electrode.TSHS
       inquire(file=TRIM(HSFile),exist=exist)
       if ( .not. exist ) then
          call die(trim(LR)//" electrode file does not exist. &
               &Please create electrode '"//trim(HSFile)//"' first.")
       end if
       ! Read in the number of atoms in the HSfile
       call ts_read_TSHS_na(HSFile,tmp_NUsedAtoms)

       if ( NUsedAtoms < 0 ) then
          NUsedAtoms = tmp_NUsedAtoms
       else if ( NUsedAtoms == 0 ) then
          if(IONode) &
               write(*,*) "You need at least one atom in the electrode."
          call die("None atoms requested for electrode calculation.")
       else if ( tmp_NUsedAtoms < NUsedAtoms ) then
          if (IONode) then
             write(*,*) "# of requested atoms is larger than available."
             write(*,*) "Requested: ",NUsedAtoms
             write(*,*) "Available: ",tmp_NUsedAtoms
          end if
          call die("Error on requested atoms.")
       end if

       ! We have determined the number of atoms in the 
       ! TSHS file
       ! Read in lasto to determine the number of orbitals 
       ! used in the electrode
       allocate(lasto(0:tmp_NUsedAtoms))
       call ts_read_TSHS_lasto(HSFile,tmp_NUsedAtoms,lasto)
       NUsedOrbs = 0
       if ( LR == 'Left' ) then
          ! We use the first atoms
          do i = 1 , NUsedAtoms
             NUsedOrbs = NUsedOrbs + lasto(i)-lasto(i-1)
          end do
       else
          ! We use the last atoms
          do i = tmp_NUsedAtoms - NUsedAtoms + 1 , tmp_NUsedAtoms
             NUsedOrbs = NUsedOrbs + lasto(i)-lasto(i-1)
          end do
       end if
    end if

  end subroutine check_HSfile
  
end subroutine read_ts_options

end module m_ts_options
