program tshs2tshs
  
  use f2kcli
  use units
  use parallel
  use precision, only : dp
  use m_hs_matrix
  use m_ts_io
  use m_ts_io_version
  use geom_helper, only : ucorb

  implicit none

  ! strings used to print out energies...
  character(len=500) :: filein, fileout, arg, tmp
  integer :: vin, vout

  integer :: iarg, narg
  logical :: exists, force

  ! ******************* TSHS ********************
  logical :: onlyS
  logical :: Gamma, TSGamma
  real(dp) :: ucell(3,3)
  integer :: na_u, no_l, no_u, no_s, maxnh, nspin,nsc(3)
  real(dp), pointer :: xa(:,:) ! (3,na_u)
  integer, pointer :: numh(:), listhptr(:) ! (no_u)
  integer, pointer :: listh(:) !(maxnh)
  real(dp), pointer :: xij(:,:) ! (3,maxnh)
  integer, pointer :: lasto(:) ! (0:na_u) 
  real(dp), pointer :: H(:,:), S(:) !(maxnh,nspin),(maxnh)
  integer :: kscell(3,3)
  real(dp) :: kdispl(3)
  real(dp) :: Ef, Qtot, Temp
  integer :: istep, ia1
  ! not really part of TSHS anymore
  integer, allocatable :: iza(:)
  ! *********************************************
  integer :: n_nsc(3), i
  integer, allocatable :: indxuo(:) ! (no_s) 
  
  force = .false.
  IONode = .true.
  Node = 0
  Nodes = 1

  ! If not gamma, then this will most likely fail...
  n_nsc = 1
  ! Default version out file to be 0 (old format)
  vout = 0

  ! Here we start the routines
  filein  = 'none'
  fileout = 'none'
  narg = command_argument_count()
  iarg = 1
  do while( iarg <= narg )
     arg = ' '
     call get_command_argument(iarg,arg)
     select case ( arg )
     case ( '-v' , '--version' )
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        read(arg,'(i10)') vout
     case ( '-f' , '--force' )
        force = .true.
     case ( '-o', '--out', '-out' )
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        fileout = arg
     case ( '-h', '--help', '-help' )
        call help
     case ( '-sA' , '-sB' , '-sC' )
        iarg = iarg + 1
        call get_command_argument(iarg,tmp)
        if (      arg(3:3) == 'A' ) then
           read(tmp,'(i10)') n_nsc(1)
        else if ( arg(3:3) == 'B' ) then
           read(tmp,'(i10)') n_nsc(2)
        else if ( arg(3:3) == 'C' ) then
           read(tmp,'(i10)') n_nsc(3)
        end if
     case default
        if ( arg(1:1) == '-' ) then
           write(0,'(a)') 'Either of the two errors has been encountered for the option "'//trim(arg)//'":'

           write(0,'(a)') ' 1) The option is not recognised'
           write(0,'(a)') ' 2) Input fdf cannot start with a hyphen "-"'
           call nl(0)
           call help
        end if
        if ( filein /= 'none' ) then
           fileout = arg
        else
           filein = arg
        end if
     end select
     iarg = iarg + 1
  end do
  if ( filein == 'none' ) then
     write(0,'(a)') 'Could not find input file on the command line'
     call nl(0)
     call help
  end if
  if ( fileout == 'none' ) then
     write(0,'(a)') 'Could not find output file on the command line'
     call nl(0)
     call help
  end if

  ! check whether the file exists
  inquire(file=filein,exist=exists)
  if (.not. exists ) then
     write(0,'(a)') 'Input file does not exist...'
     call nl(0)
     call help
  end if
  inquire(file=fileout,exist=exists)
  if ( exists ) then
     write(0,'(a)') 'Out put file already exist... You cannot overwrite...'
     call nl(0)
     call help
  end if

  ! get the input version
  vin = TSHS_version(filein)

  select case ( vin )
  case ( 0, 1 )
     ! Do nothing
  case default
     write(0,'(a)') 'Version not recognized by this program, please update...'
     write(0,'(a)') 'Versions allowed are: 0,1'
     write(0,'(a,i0,a)') 'Found version: ',vin,' in file: '//trim(filein)
     call nl(0)
     call help
  end select

  select case ( vout )
  case ( 0, 1 )
     ! Do nothing
  case default
     write(0,'(a)') 'Version not recognized by this program, please update...'
     write(0,'(a)') 'Versions allowed are: 0,1'
     write(0,'(a,i0,a)') 'Found version: ',vin,' in file: '//trim(fileout)
     call nl(0)
     call help
  end select

  if ( vin == vout .and. .not. force ) then
     write(0,'(a)') 'There is no need to convert the file.'
     write(0,'(a)') 'The version is the same...'
     call nl(0)
     call help
  end if

  write(*,'(a)') 'Reading in '//trim(filein)

  ! Read in TSHS
  call ts_read_TSHS(filein,onlyS,Gamma,TSGamma, &
       ucell, nsc, na_u, no_l, no_u, no_s, maxnh, nspin, &
       kscell, kdispl, &
       xa, lasto, &
       numh, listhptr, listh, xij , &
       H, S, Ef, Qtot, Temp, istep, ia1)
  if ( all(nsc == 0) ) then
     nsc = n_nsc
  end if

  allocate(indxuo(no_s))
  do i = 1 , no_s
     indxuo(i) = ucorb(i,no_u)
  end do
  write(*,'(a)') 'Writing to '//trim(fileout)

  if ( vout == 0 ) then
     allocate(iza(na_u))
     iza = 0
     call write_TSHS_0(fileout, &
          onlyS, Gamma, TSGamma, &
          ucell, na_u, no_l, no_u, no_s, maxnh, nspin,  &
          kscell, kdispl, &
          xa, iza, lasto, &
          numh, listhptr, listh, xij, &
          H, S, Ef, Qtot, Temp, istep, ia1)
     deallocate(iza)
  else if ( vout == 1 ) then
     ! Hopefully nsc is correct
     write(*,'(a)') 'If this fails, try forcing the supercell with -s[ABC]'
     i = no_s / no_u
     write(*,'(a,i0)') 'Hint, total supercells: ',i
     if ( product(nsc) /= i ) then
        write(*,'(a)') 'You are trying a different number &
             &of supercells, please do not do that...'
        write(*,'(a,3(tr1,i0),a,i0)') 'nsc(3) = ',nsc,'. prod = ',product(nsc)
        stop 'Stopping!'
     end if
     call ts_write_tshs(fileout, onlyS, Gamma, TSGamma, &
          ucell, nsc, na_u, no_l, no_u, no_s, maxnh, nspin, &
          kscell, kdispl, &
          xa, lasto, &
          numh, listhptr, listh, xij, indxuo, H, S, Ef, &
          Qtot, Temp, istep, ia1)
  end if

  deallocate(indxuo)
  deallocate(xa,numh,listhptr,listh)
  deallocate(xij,lasto,H,S)
  write(*,'(a)') 'File written.'

contains
  
  subroutine nl(u)
    integer, intent(in), optional :: u
    if ( present(u) ) then
       write(u,*)
    else
       write(*,*)
    end if
  end subroutine nl

  subroutine help()
    write(0,'(a)') 'Helps converting an old TranSIESTA input to the new format'
    write(0,'(a)') 'Options:'
    write(0,'(a)') '  -v|--version <integer>:'
    write(0,'(a)') '          defines the output version of the TSHS file:'
    write(0,'(a)') '            -v 0 is the old TSHS format'
    write(0,'(a)') '            -v >0 is a newer TSHS format'
    write(0,'(a)') '  -sA|-sB|-sC <integer> > 0:'
    write(0,'(a)') '          defines the number of supercell repetitions of the TSHS file:'
    write(0,'(a)') '            -sA 1 means a supercell repetition 0 in the first lattice vector'
    write(0,'(a)') '               You can read of the supercell of in the output of SIESTA'
    write(0,'(a)') '  -o|--out <filename>:'
    write(0,'(a)') '          The output file name (must not exist)'
    write(0,'(a)') ' -h|--help : this help'
    stop
  end subroutine help

end program tshs2tshs    

