! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_tbt_options
! ************************
! * SIESTA modules       *
! ************************
  use precision, only : dp
  use diagmemory, only : MemoryFactor

  implicit none

  PUBLIC
  SAVE

! ##########################
! # SIESTA options...      #
! ##########################
  character(len=200) :: sname

! ################################################
! #                                              #
! #  Similar options to those found in           #
! #  TranSIESTA. For the majority of the         #
! #  options TS.<option> there is an equivalent  #
! #  TS.TBT.<option>. The TBT takes precendence  #
! #  if both are found.                          #
! #                                              #
! ################################################
  
  logical  :: UseBulk      ! Use Bulk Hamiltonian in Electrodes
  real(dp) :: VoltFDF      ! Bias applied, Internally Volt=voltfdf/eV. 
  real(dp) :: VoltL        ! Bias on the left electrode   (  .5 * VoltFDF )
  real(dp) :: VoltR        ! Bias on the right electrode  ( -.5 * VoltFDF )
  logical  :: IsVolt       ! Has the value VoltFDF > 0.001/eV
  real(dp) :: kT           ! Electronic temperature
  real(dp) :: GFEta        ! Imaginary part of the Bias Contour  
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
  logical :: ElecValenceBandBot ! Calculate Electrode valence band bottom when creating electrode GF
  logical :: ReUseGF        ! Calculate the electrodes GF

! ################################################
! #                                              #
! #  Default values to those found in            #
! #  TranSIESTA.                                 #
! #                                              #
! ################################################
  logical, parameter :: UseBulk_def = .true.
  real(dp), parameter :: VoltFDF_def = 0._dp      ! in Ry
  real(dp), parameter :: kT_def = 1.9e-3_dp       ! in Ry
  real(dp), parameter :: GFEta_def = 0.000001_dp  ! in Ry
  integer, parameter :: NBufAt_def = 0
  integer, parameter :: NRepA_def = 1
  integer, parameter :: NUsedAtoms_def = -1
  character(33), parameter :: GFTitle_def = 'Generated GF file'
  character(33), parameter :: HSFile_def = 'NOT REQUESTED'
  logical, parameter :: ElecValenceBandBot_def = .true.
  logical, parameter :: ReUseGF_def = .true.

! ################################################
! #                                              #
! #  Specific options for TBTrans                #
! #                                              #
! ################################################
  character(200) :: HSFile ! The scattering region TSHS file
  real(dp)       :: Emin ! The minimum energy to evaluate the transmission on
  real(dp)       :: Emax ! The maximum energy to evaluate the transmission on
!                          Thus the evaluation will be in the range [Emin; Emax]
  integer        :: NPoints  ! the number of devisions on the energy range
  integer        :: NeigCh   ! No. eigenchannels to calculate
  logical        :: CalcIeig ! Calculate the eigenvalues in the projected isolated region
  integer        :: IsoAt1   ! The first atom in the isolated region
  integer        :: IsoAt2   ! The last atom in the isolated region
  logical        :: SpinPol  ! Spin polarized calculation?
  logical        :: CalcCOOP ! Do the COOP curves
  logical        :: AlignScat ! Align the scattering region with the left electrode (only by the first onsite element)
  logical        :: CalcAtomPDOS ! Calculate the DOS on the projected atoms
  logical        :: RemUCellDistances ! Remove the Ucell distances when calculating the transmission...

! ################################################
! #                                              #
! #  Default options for TBTrans                 #
! #                                              #
! ################################################
  real(dp), parameter :: Emin_def = -.2_dp
  real(dp), parameter :: Emax_def =  .2_dp
  integer, parameter  :: NPoints_def = 100
  integer, parameter  :: NeigCh_def = 5
  logical, parameter  :: SpinPol_def = .false.
  logical, parameter  :: CalcIeig_def = .false.
  logical, parameter  :: CalcCOOP_def = .false.
  logical, parameter  :: AlignScat_def = .false.
  logical, parameter  :: CalcAtomPDOS_def = .false.
  logical, parameter  :: RemUCellDistances_def = .false.

CONTAINS
  
  ! Read in the tbtrans options
  subroutine read_tbt_options()

! ************************
! * SIESTA modules       *
! ************************
    use fdf,            only : leqi, fdf_defined, fdf_deprecated, fdf_obsolete
    use parallel,       only : IOnode, Nodes, operator(.PARCOUNT.)
    use m_fdf_global,   only : fdf_global_get
    use units,          only : eV
    use files,          only : slabel
    use sys,            only : die
    use m_ts_io       , only : ts_read_TSHS_na
    use m_ts_io       , only : ts_read_TSHS_lasto
#ifdef MPI
    use mpi_siesta, only: MPI_Bcast, MPI_character, MPI_Comm_World
#endif


! Internal Variables
    character(len=200) :: chars
    logical :: exist ! Check file existance for files requested
    character(len=200) :: paste
    integer :: na_u, tmp, i
    external :: paste
#ifdef MPI
    integer :: MPIerror
#endif
    
    if (IOnode) then
       write(*,*)
       write(*,'(2a)') 'tbt_read_options: ', repeat('*', 62)
    end if

    ! Show the deprecated and obsolete labels
    call fdf_deprecated('TS.TBT.DoCOOP','TS.TBT.COOP')
    call fdf_deprecated('TS.CalcGF','TS.TBT.ReUseGF')

    
! Reading from fdf ... This is needed for using 'cdiag'
    call fdf_global_get(MemoryFactor,'Diag.Memory', 1.0_dp )
    
    call fdf_global_get(Emin,'TS.TBT.Emin',Emin_def,'Ry')
    call fdf_global_get(Emax,'TS.TBT.Emax',Emax_def,'Ry')
    call fdf_global_get(GFEta,'TS.TBT.Eta',GFeta_def,'Ry')
    call fdf_global_get(VoltFDF,'TS.Voltage',VoltFDF_def,'Ry') 
    IsVolt = dabs(VoltFDF) > 0.001_dp/eV
    ! Assign the fermi shifts in the electrodes
    VoltL =  .5_dp*VoltFDF
    VoltR = -.5_dp*VoltFDF
    call fdf_global_get(NPoints,'TS.TBT.NPoints',NPoints_def)
    call fdf_global_get(Neigch,'TS.TBT.NEigen',Neigch_def)
    call fdf_global_get(SpinPol,'SpinPolarized',SpinPol_def)
    call fdf_global_get(UseBulk,'TS.UseBulkInElectrodes',UseBulk_def)
    call fdf_global_get(NBufAtL,'TS.BufferAtomsLeft',NBufAt_def)
    call fdf_global_get(NBufAtR,'TS.BufferAtomsRight',NBufAt_def)
    call fdf_global_get(kT,'ElectronicTemperature',kT_def,'Ry')
    call fdf_global_get(GFTitle,'TS.TBT.GFTitle',GFTitle_def)
    chars = trim(slabel)//'.TBTGFL'
    call fdf_global_get(GFFileL,'TS.TBT.GFFileLeft',trim(chars))
    chars = trim(slabel)//'.TBTGFR'
    call fdf_global_get(GFFileR,'TS.TBT.GFFileRight',trim(chars))
    call fdf_global_get(ReUseGF,'TS.TBT.ReUseGF',ReUseGF_def)
    ! This needs a two way entrance (in TranSIESTA it really doesn't matter.
    ! In TBTrans it can be used to check for Emin against the valence band bottom
    call fdf_global_get(ElecValenceBandBot,'TS.TBT.CalcElectrodeValenceBandBottom', &
         ElecValenceBandBot_def)
    chars = paste(slabel,'.TSHS')
    call fdf_global_get(HSFile,'TS.TBT.HSFile',chars)
    ! Check for file existance
    inquire(file=TRIM(HSFile),exist=exist)
    if ( .not. exist ) then
       call die("Scattering region does not exist. &
            &Please create scattering region file '"//TRIM(HSFile)//"' first.")
    end if
    ! Read in total number of atoms in the TSHS file!
    call ts_read_TSHS_na(HSFile,na_u)

    ! Read electrode options
    call fdf_global_get(HSFileL,'TS.HSFileLeft',HSFile_def)
    call fdf_global_get(NUsedAtomsL,'TS.NumUsedAtomsLeft',NUsedAtoms_def)
    call check_HSfile('Left',HSFileL,NUsedAtomsL,NUsedOrbsL)
    call fdf_global_get(NRepA1L,'TS.ReplicateA1Left',NRepA_def)
    call fdf_global_get(NRepA2L,'TS.ReplicateA2Left',NRepA_def)
    if ( NRepA1L < 1 .or. NRepA2L < 1 ) &
         call die("Repetition in left electrode must be >= 1.")

    call fdf_global_get(HSFileR,'TS.HSFileRight',HSFile_def)
    call fdf_global_get(NUsedAtomsR,'TS.NumUsedAtomsRight',NUsedAtoms_def)
    call check_HSfile('Right',HSFileR,NUsedAtomsR,NUsedOrbsR)
    call fdf_global_get(NRepA1R,'TS.ReplicateA1Right',NRepA_def)
    call fdf_global_get(NRepA2R,'TS.ReplicateA2Right',NRepA_def)
    if ( NRepA1R < 1 .or. NRepA2R < 1 ) &
         call die("Repetition in right electrode must be >= 1.")

    call fdf_global_get(CalcIeig,'TS.TBT.CalcIeig',CalcIeig_def)

    tmp = 1+NBufAtL+NUsedAtomsL*NRepA1L*NRepA2L
    call fdf_global_get(IsoAt1,'TS.TBT.PDOSFrom',tmp)
    if ( IsoAt1 < tmp ) then
       call die("Requested PDOS 1 atom is outside of contact region. &
            &Please choose an atom within the device.")
    end if
    tmp = na_u-NBufAtR-NUsedAtomsR*NRepA1R*NRepA2R
    call fdf_global_get(IsoAt2,'TS.TBT.PDOSTo',tmp)
    if ( IsoAt2 > tmp ) then
       call die("Requested PDOS 2 atom is outside of contact region. &
            &Please choose an atom within the device.")
    end if

    call fdf_global_get(CalcCOOP,'TS.TBT.COOP',CalcCOOP_def)
    call fdf_global_get(AlignScat,'TS.TBT.AlignOnSite',AlignScat_def)
    call fdf_global_get(CalcAtomPDOS,'TS.TBT.AtomPDOS',CalcAtomPDOS_def)
    call fdf_global_get(RemUCellDistances,'TS.TBT.RemoveUnitCellDistance',RemUCellDistances_def)

! Output Used Options in OUT file ....
    if (ionode) then
       write(*,6) 'TBTrans Voltage                               =', voltfdf/eV,' V'
       write(*,6) 'TBTrans Emin                                  =', Emin/eV,' eV'
       write(*,6) 'TBTrans Emax                                  =', Emax/eV,' eV'
       write(*,1) 'Bulk Values in Electrodes                     =', UseBulk
       write(*,5) 'Buffer Atoms in Left electrode                =', NBufAtL
       write(*,5) 'Buffer Atoms in Right electrode               =', NBufAtR
       write(*,5) 'Points on the energy contour                  =', NPoints
       write(*,7) 'GFEta                                         =', GFEta,' Ry'
       write(*,6) 'Electronic Temperature                        =', kT, ' Ry'
       write(*,1) 'Calculate band bottom in elecrodes            =', ElecValenceBandBot
       write(*,10)'GF title                                      =', TRIM(GFTitle)
       write(*,10)'Left GF File                                  =', TRIM(GFFileL)
       write(*,10)'Right GF File                                 =', TRIM(GFFileR)
       write(*,1) 'Re-use the GF files if they exists            =', ReUseGF
       write(*,10)'Scattering region TSHS file                   =', TRIM(HSFile)
       write(*,10)'Left electrode TSHS file                      =', TRIM(HSFileL)
       write(*,5) '# atoms used in left elec.                    = ', NUsedAtomsL
       write(*,'(a,i3,'' X '',i3)') &
                  'Left elec. repetition A1/A2                   = ', NRepA1L,NRepA2L
! Check existance for right Electrode.TSHS
       write(*,10)'Right electrode TSHS file                     =', TRIM(HSFileL)
       write(*,5) '# atoms used in right elec.                   = ', NUsedAtomsL
       write(*,'(a,i3,'' X '',i3)') &
                  'Right elec. repetition A1/A2                  = ', NRepA1L,NRepA2L
       write(*,'(a,''['',i5,'';'',i5,'']'')') &
                  'Projected region                              = ', IsoAt1,IsoAt2
       write(*,1) 'Calculate DOS on projected atoms              = ',CalcAtomPDOS
       write(*,1) 'Calculate COOP                                = ',CalcCOOP
       write(*,1) 'Align the Hamiltonian with the electrode      = ',AlignScat
       write(*,1) 'Remove inner-cell distances in the Hamiltonian= ',RemUCellDistances
       if ( AlignScat ) then
          call die("TBtrans is currently not implented to align the scattering &
               &region and the electrodes.")
       end if
    end if
 
    if (IOnode) then
       write(*,'(2a,/)') 'tbt_read_options: ', repeat('*', 62)
    end if

! Print out message if the number of contour points are not 
! divisable by the number of Nodes
    if ( IONode .and. mod(NPoints,Nodes) /= 0 ) then
       write(*,*) "NOTICE: Transport energy points are not"
       write(*,*) "        divisable by the number of nodes."
       write(*,*) "        Better scalability is achived by changing:"
       write(*,*) "          - TS.TBT.NPoints"

! Calculate optimal number of energy points
       i = NPoints
       write(*,'(t10,a,i4)') "Used # of energy points   : ",i
       i = Nodes .PARCOUNT. i
       write(*,'(t10,a,i4,tr1,a4,i3,/)') &
            "Optimal equilibrium # of energy points: ",i, &
            achar(177)//" i*",Nodes
    end if

1   format(a,4x,l1)
5   format(a,i5,a)
6   format(a,f10.4,a)
7   format(a,f12.6,a)
10  format(a,4x,a)

  contains

    subroutine check_HSfile(LR,HSFile,NUsedAtoms,NUsedOrbs)
      character(len=*), intent(in) :: LR
      character(len=*), intent(in) :: HSFile
      integer, intent(inout) :: NUsedAtoms, NUsedOrbs
      integer :: tmp_NUsedAtoms
      integer, allocatable, dimension(:) :: lasto
      logical :: exist
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

    end subroutine check_HSfile

  end subroutine read_tbt_options

end module m_tbt_options
