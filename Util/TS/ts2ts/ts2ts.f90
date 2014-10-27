! Utility to convert an FDF input in old format
! to an FDF input in new format.

program ts2ts
  
  use f2kcli
  use fdf
  use units
  use precision, only : dp
  
  implicit none

  ! strings used to print out energies...
  integer, parameter :: N_char = 50
  character(len=N_char) :: c_CCEmin, c_GFEta, c_Volt, c_dVolt
  
  real(dp) :: CCEmin, GFEta, Volt
  integer :: Nline, Ncircle, Npol, Nvolt

  logical :: IsVolt, Bulk, UpdateDMCR, ReUse

  integer :: NBufL, NBufR
  integer :: na_usedL, NRepA1L, NRepA2L
  integer :: na_usedR, NRepA1R, NRepA2R
  character(len=500) :: L_TSHS, R_TSHS, filein, arg

  integer :: iarg, narg
  logical :: exists, def_nEq


  ! Set default values (then options can overwrite them)
  Nline = -1
  Npol  = -1
  Ncircle = -1
  Nvolt = -1

  def_nEq = .false.
  ! Here we start the routines
  filein = 'none'
  narg = command_argument_count()
  iarg = 1
  do while( iarg <= narg )
     arg = ' '
     call get_command_argument(iarg,arg)
     select case ( arg )
     case ( '-neq' )
        def_nEq = .true.
     case ( '-N-pole' )
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        read(arg,'(i10)') Npol
        write(0,'(a)') '# Overwriting number of pole points'
     case ( '-N-line' )
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        read(arg,'(i10)') Nline
        write(0,'(a)') '# Overwriting number of points on Fermi-line'
     case ( '-N-circle' )
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        read(arg,'(i10)') Ncircle
        write(0,'(a)') '# Overwriting number of points on circle'
     case ( '-h', '--help', '-help' )
        call help
     case default
        if ( arg(1:1) == '-' ) then
           write(0,'(a)') 'Either of the two errors has been encountered for the option "'//trim(arg)//'":'

           write(0,'(a)') ' 1) The option is not recognised'
           write(0,'(a)') ' 2) Input fdf cannot start with a hyphen "-"'
           call nl(0)
           call help
        end if
        filein = arg
     end select
     iarg = iarg + 1
  end do
  if ( leqi(filein,'none') ) then
     write(0,'(a)') 'Could not find input file on the command line'
     call nl(0)
     call help
  end if

  ! check whether the file exists
  inquire(file=filein,exist=exists)
  if (.not. exists ) then
     write(0,'(a)') 'Input file can not be piped into this program, &
          &please supply FDF file on command line...'
     call nl(0)
     call help
  end if

  ! Initialize the fdf
  call fdf_init(filein,"ts2ts.log")

  ! Buffer atoms
  NBufL    = fdf_get('TS.BufferAtomsLeft',0)
  NBufR    = fdf_get('TS.BufferAtomsRight',0)

  ! Read in electrode stuff
  NRepA1L  = fdf_get('TS.ReplicateA1Left',1)
  NRepA2L  = fdf_get('TS.ReplicateA2Left',1)
  na_usedL = fdf_get('TS.NumUsedAtomsLeft',0)
  L_TSHS   = fdf_get('TS.HSFileLeft','none')

  NRepA1R  = fdf_get('TS.ReplicateA1Right',1)
  NRepA2R  = fdf_get('TS.ReplicateA2Right',1)
  na_usedR = fdf_get('TS.NumUsedAtomsRight',0)
  R_TSHS   = fdf_get('TS.HSFileRight','none')

  
  Volt    = fdf_get('TS.Voltage',0._dp,'Ry')
  IsVolt = abs(Volt/eV) > 0.00001_dp
  call e2a(Volt,c_Volt,force_eV=.true.)

  if ( Nline < 0 ) &
       Nline   = fdf_get('TS.ComplexContour.NLine',10)
  call assert(Nline > 0, 'Nline has to be larger than 0')
  if ( Npol < 0 ) &
       Npol    = fdf_get('TS.ComplexContour.NPoles',6)
  call assert(Npol > 0, 'NPoles has to be larger than 0')
  if ( Ncircle < 0 ) &
       Ncircle = fdf_get('TS.ComplexContour.NCircle',30)
  call assert(NCircle > 0, 'NCircle has to be larger than 0')
  if ( Nvolt < 0 ) &
       Nvolt   = fdf_get('TS.biasContour.NumPoints',5)
  if ( IsVolt ) then
     call assert(Nvolt > 0, 'NumPoints has to be larger than 0')
     call e2a(abs(Volt)/Nvolt,c_dVolt,force_eV=.true.)
     if ( def_nEq ) then
        c_dVolt = '0.01 eV'
     end if
  else
     ! for no bias, we can only guess what the user would do
     ! Hence we default it to the correct thing
     c_dVolt = '0.01 eV'
  end if
  ! Read in contour stuff...
  CCEmin  = fdf_get('TS.ComplexContourEmin',-3._dp,'Ry')
  call e2a(CCEmin,c_CCEmin)
  GFEta   = fdf_get('TS.biasContour.Eta',0.000001_dp,'Ry')
  call e2a(GFEta,c_GFEta,prec=10,force_eV=.false.)

  Bulk       = fdf_get('TS.UseBulkInElectrodes',.true.)
  ! in the original implementation this meant only update the central region
  ! as such I thought the key-word could be interpreted as DMCR == DM-crossterms
  ! So I have "spelled it out"
  UpdateDMCR = .not. fdf_get('TS.UpdateDMCROnly',.true.)
  ReUse      = .not. fdf_get('TS.calcGF',.true.)
  ReUse      = fdf_get('TS.ReUseGF',ReUse)

  ! Start to write out the information...
  write(*,'(a,a)')'TS.Voltage ',c_Volt
  call nl

  ! Start with writing out the buffer atoms:
  if ( NBufL + NBufR > 0 ) then
     ! Write block
     call sblock('TS.Atoms.Buffer')
     if ( NBufL > 1 ) then
        write(*,'(tr2,a,i0,a,i0)') 'atom from ',1,' to ',NBufL
     else if ( NBufL == 1 ) then
        write(*,'(tr2,a,i0,a,i0)') 'atom 1'
     end if
     if ( NBufR > 1 ) then
        write(*,'(tr2,a,i0,a,i0)') 'atom from ',-NBufR,' to ',-1
     else if ( NBufR == 1 ) then
        write(*,'(tr2,a,i0,a,i0)') 'atom -1'
     end if
     call eblock('TS.Atoms.Buffer')
     call nl()
  end if

  ! Write out the chemical potentials, they are pretty standard... :)
  call sblock('TS.ChemPots')
  write(*,'(tr2,a)') 'Left'
  write(*,'(tr2,a)') 'Right'
  call eblock('TS.ChemPots')
  call nl
  call wchem('Left',.true.)
  call wchem('Right',.false.)

  call nl()
  ! Write out the generic options for the electrodes
  if ( Bulk ) then
     write(*,'(a)') 'TS.Elecs.Bulk true'
  else
     write(*,'(a)') 'TS.Elecs.Bulk false'
  end if
  if ( UpdateDMCR ) then
     write(*,'(a)') 'TS.Elecs.DM.Update cross-terms'
  else
     write(*,'(a)') 'TS.Elecs.DM.Update none'
  end if
  if ( ReUse ) then
     write(*,'(a)') 'TS.Elecs.GF.ReUse true'
  else
     write(*,'(a)') 'TS.Elecs.GF.ReUse false'
  end if

  call sblock('TS.Elecs')
  write(*,'(tr2,a)') 'Left'
  write(*,'(tr2,a)') 'Right'
  call eblock('TS.Elecs')
  call nl
  ! Write out the electrodes...
  call welec('Left',.true.  , L_TSHS,NBufL,NRepA1L,NRepA2L,na_usedL)
  call nl
  call welec('Right',.false., R_TSHS,NBufR,NRepA1R,NRepA2R,na_usedR)
  call nl

  ! Start by writing the contours
  write(*,'(a,i0)') 'TS.Contours.Eq.Pole.N ',Npol
  call wcont('c-Left','circle','g-legendre', &
       trim(c_CCEmin)//' + V/2','-10. kT + V/2',Ncircle)
  call wcont('t-Left','tail','g-fermi', 'prev','inf',Nline)
  call wcont('c-Right','circle','g-legendre', &
       trim(c_CCEmin)//' - V/2','-10. kT - V/2',Ncircle)
  call wcont('t-Right','tail','g-fermi', 'prev','inf',Nline)
  call nl

  ! Now we need to write out the non-equilibrium contour...
  write(*,'(a,a)')'TS.Contours.nEq.Eta ',c_GFEta
  call sblock('TS.Contours.Bias.Window')
  write(*,'(tr2,a)') 'neq'
  call eblock('TS.Contours.Bias.Window')

  ! Now we better tell the user that delta 0.01 eV is
  ! the best thing!
  if ( (.not. def_nEq ).and. abs(Volt) / Nvolt > 0.05_dp * eV ) then
     write(0,'(a)') '# Selected non-equilibrium integration'
     write(0,'(a)') '# integration splitting is larger than 0.05 eV.'
     write(0,'(a)') '# New default setting is: 0.01 eV.'
     write(0,'(a)') '# Consider changing the non-equilibrium "delta <E>" to 0.01 eV.'
     write(*,'(a)') '# Selected non-equilibrium integration'
     write(*,'(a)') '# integration splitting is larger than 0.05 eV.'
     write(*,'(a)') '# New default setting is: 0.01 eV.'
     write(*,'(a)') '# Consider changing the below "delta <E>" to 0.01 eV.'

  end if
  call wcont('Bias.Window.neq','line','simpson-mix', &
       '-|V|/2','|V|/2',cD=c_dVolt)
  call nl
  call sblock('TS.Contours.Bias.Tail')
  write(*,'(tr2,a)') 'neq-tail'
  call eblock('TS.Contours.Bias.tail')
  call wcont('Bias.Tail.neq-tail','tail','simpson-mix', &
       '0. kT','12. kT',cD=c_dVolt)

contains
  
  subroutine e2a(e,a,force_eV,prec)
    real(dp), intent(in) :: e
    character(len=N_char),intent(out) :: a
    integer, intent(in), optional :: prec

    logical, intent(in), optional :: force_eV
    character(len=N_char) :: ctmp
    integer :: N_Ry, N_eV, lprec
    real(dp), parameter :: cr = 0.001_dp
    real(dp) :: tmp

    lprec = 5
    if ( present(prec) ) lprec = prec

    write(ctmp,'(f30.20)') e - nint(e)
    N_Ry = ccount0(ctmp)
    write(ctmp,'(f30.20)') e / eV - nint(e/eV)
    N_eV = ccount0(ctmp)

    a = ' '
    if ( present(force_eV) ) then
       if ( force_eV ) then
          N_eV = N_Ry + 1
       else
          N_Ry = N_eV + 1
       end if
    end if
    write(ctmp,'(a,i0,a,i0,a)')'(f',lprec+5,'.',lprec,',tr1,a)'
    if ( N_eV > N_Ry ) then
       ! correct
       if ( abs(e/eV - nint(e/eV)) < cr .and. abs(e/eV) > cr ) then
          tmp = nint(e / eV)
       else
          tmp =      e / eV
       end if
       ! the unit must be given in eV
       write(a,ctmp) tmp,'eV'
    else
       if ( abs(e - nint(e)) < cr .and. abs(e) > cr ) then
          tmp = nint(e)
       else
          tmp =      e
       end if
       write(a,ctmp) tmp,'Ry'
    end if
    if ( ccount0(a) > lprec ) then
       write(0,'(a)')'###'
       write(0,'(2a)')'### Please check your input, a number might be interpreted as 0: ',trim(a)
       write(a,'(g20.10,a)') e,' Ry'
       write(0,'(2a)')'### Will revert to this: ',trim(a)

       write(0,'(a)')'###'
    end if
  end subroutine e2a

  function ccount0(str) result(N)
    character(len=*), intent(in) :: str
    integer :: i, N, N0
    N  = 0
    N0 = 0
    do i = 1 , len(str)
       if ( str(i:i) == '0' .and. N == 0 ) then
          N0 = N0 + 1
       else if ( str(i:i) == '-' .or. &
            str(i:i) == ' ' .or. &
            str(i:i) == '.' ) then
          ! do nothing
       else if ( str(i:i) == '9' .and. N0 == 1 ) then
          N = N + 1
       else
          exit
       end if
    end do
    N = N + N0
  end function ccount0
  
  subroutine nl(u)
    integer, intent(in), optional :: u
    if ( present(u) ) then
       write(u,*)
    else
       write(*,*)
    end if
  end subroutine nl
  subroutine sblock(block)
    character(len=*), intent(in) :: block
    write(*,'(2a)') '%block ',trim(block)
  end subroutine sblock
  subroutine eblock(block)
    character(len=*), intent(in) :: block
    write(*,'(2a)') '%endblock ',trim(block)
  end subroutine eblock
  
  subroutine wchem(name,is_Left)
    character(len=*), intent(in) :: name
    logical, intent(in) :: is_Left
    call sblock('TS.ChemPot.'//trim(name))
    if ( is_Left ) then
       write(*,'(tr2,a)') 'mu V/2'
    else
       write(*,'(tr2,a)') 'mu -V/2'
    end if
    write(*,'(tr2,a)') 'contour.eq'
    write(*,'(tr4,a)') 'begin'
    write(*,'(tr6,a)') 'c-'//trim(name)
    write(*,'(tr6,a)') 't-'//trim(name)
    write(*,'(tr4,a)') 'end'
    call eblock('TS.ChemPot.'//trim(name))
  end subroutine wchem

  subroutine welec(name,is_Left,TSHS,NBuf,RepA1,RepA2,na)
    character(len=*), intent(in) :: name,TSHS
    logical, intent(in) :: is_Left
    integer, intent(in) :: NBuf, RepA1, RepA2, na

    call sblock('TS.Elec.'//trim(name))
    write(*,'(tr2,2a)') 'TSHS ',trim(TSHS)
    if ( is_Left ) then
       write(*,'(tr2,a)') 'chem-pot Left'
       write(*,'(tr2,a)') 'semi-inf-dir -a3'
       write(*,'(tr2,a,i0)') 'elec-pos begin ',NBuf+1
    else
       write(*,'(tr2,a)') 'chem-pot Right'
       write(*,'(tr2,a)') 'semi-inf-dir +a3'
       write(*,'(tr2,a,i0)') 'elec-pos end ',-NBuf-1
    end if
    if ( na > 0 ) then
       write(*,'(tr2,a,i0)') 'used-atoms ', na
    end if
    ! Only print out repetitions if larger than one
    if ( RepA1 > 1 ) then
       write(*,'(tr2,a,i0)') 'rep-a1 ', RepA1
    end if
    if ( RepA2 > 1 ) then
       write(*,'(tr2,a,i0)') 'rep-a2 ', RepA2
    end if
    call eblock('TS.Elec.'//trim(name))
  end subroutine welec

  subroutine wcont(name,part,method,cEmin, cEmax,N,cD)
    character(len=*), intent(in) :: name, part, method, cEmin, cEmax
    integer, intent(in), optional :: N
    character(len=*), intent(in), optional :: cD
    call sblock('TS.Contour.'//trim(name))
    write(*,'(tr2,2a)')'part ',trim(part)
    write(*,'(tr2,4(tr1,a))')'from',trim(cEmin),'to',trim(cEmax)
    if ( present(N) )then
       write(*,'(tr4,a,i0)')'points ',N
    else
       write(*,'(tr4,2a)')'delta ',trim(cD)
    end if
    write(*,'(tr4,2(tr1,a))')'method',trim(method)
    call eblock('TS.Contour.'//trim(name))
  end subroutine wcont

  subroutine assert(bool,msg)
    logical, intent(in) :: bool
    character(len=*), intent(in) :: msg
    if ( .not. bool ) then
       write(0,'(a)')'### ERROR ###'
       write(0,'(2a)')'### ',trim(msg)
       write(0,'(a)')'### ERROR ###'
       stop
    end if
  end subroutine assert

  subroutine help()
    write(0,'(a)') 'Helps converting an old TranSIESTA input to the new format'
    write(0,'(a)') 'Options:'
    write(0,'(a)') '  -neq      : overrides the number of non-equilibrium contour'
    write(0,'(a)') '              points to be "delta 0.01 eV" so it is independently'
    write(0,'(a)') '              set for all bias voltages.'
    write(0,'(a)') '  -N-pole <int> | -N-circle <int> | -N-line <int> :'
    write(0,'(a)') '              overrides the number of points on the equivalent contour'
    write(0,'(a)') ' <fdf>      : the input fdf file that needs conversion.'
    write(0,'(a)') ' -h         : this help.'
    stop
  end subroutine help

end program ts2ts
    

