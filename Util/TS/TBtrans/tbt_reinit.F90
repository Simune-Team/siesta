! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
! Transferred to be used in the TBTrans utility.
!
subroutine tbt_reinit( sname , slabel ) 

! Subroutine to initialise the reading of the data for SIESTA 
!
!     It uses the FDF (Flexible Data Format) package 
!     of J.M.Soler and A.Garcia
!
! Taken from redata. Writen by P.Ordejon, December'96
! **************************** OUTPUT *********************************
! character    slabel      : System Label (to name output files)
! character(len=*) sname       : System Name
! **********************************************************************

!
!  Modules
!
  use parallel,    only : Node
  use fdf
  use m_verbosity

  implicit none

  character(len=*), intent(out) :: sname, slabel

!  Internal variables .................................................
  character(len=50) :: filein, fileout, string

  integer ::  count, length, lun, lun_tmp, iostat
  character(len=256) :: line

  logical :: debug_input, file_exists

! Print Welcome and Presentation .......................................
!     Non-master mpi-processes receive a copy of all the
!     pre-processed fdf input information (recursively
!     including files when necessary), and dump it on a
!     text file with name "fdf_input.<ProcessNumber>".
!     They then read from this file to construct a memory
!     image of the information.
!
!     The master process creates its memory image directly
!     from the standard fdf file (renamed as INPUT_TMP.$$,
!     as usually done for historical reasons (see below)).
!
  filein = "fdf_input" 

  if (Node.eq.0) then
     write(6,'(/a)') &
          '                           ************************ '
#ifdef TBT_PHONON
     write(6,'(a)') &
          '                           *  WELCOME TO PHtrans  * '
#else
     write(6,'(a)') &
          '                           *  WELCOME TO TBtrans  * '
#endif
     write(6,'(a)') &
          '                           ************************ '
! ..................
!
!       Set name of file to read from. Done only
!       in the master node.
!
!
!     Choose proper file for fdf processing
!     (INPUT_DEBUG if it exists or "standard input",
!      processed and dumped to a temporary file)
!
     inquire(file='INPUT_DEBUG',exist=debug_input)
     if (debug_input) then
        write(*,'(a)') 'WARNING: ' // &
             'TBTrans is reading its input from file INPUT_DEBUG'
        filein = 'INPUT_DEBUG'

#ifndef NO_F2003
     else if ( command_argument_count() > 0 ) then

        ! Get file-name from input line
        filein = ' '
        call get_command_argument(1,filein,length)
        if ( length > len(filein) ) then
           call die('The argument is too long to be retrieved, please &
                &limit to 50 characters for the input file')
        end if
        inquire(file=filein,exist=debug_input)
        if ( .not. debug_input ) then
           call die('Input file '//trim(filein)//' does not exist? Have &
                &you specified the wrong file-name?')
        end if

        write(*,'(/,2a)') 'reinit: Reading from ',trim(filein)


#endif
     else
!
!          Read from standard input (dumped to a temp file)
!
        write(*,'(/a)') 'reinit: Reading from standard input'
        lun = 5
        call io_assign(lun_tmp)
        do  ! make sure we get a new file
           call system_clock( count )
           write(string,*) count
           filein = 'INPUT_TMP.'//adjustl(string)
           inquire( file=filein, exist=file_exists )
           if (.not.file_exists) exit
        end do
!
        open(lun_tmp,file=filein, &
             form='formatted',status='replace')
        rewind(lun_tmp)
        write(*,"(a,23('*'),a,28('*'))") &
             '***', ' Dump of input data file '
!
        do
           read(lun,iostat=iostat,fmt='(a)') line
           if (iostat /= 0 ) exit
           length = len_trim(line)
           if (length /= 0) then
              write(*,'(a)') line(1:length)
              if (.not. debug_input) then
                 write(lun_tmp,'(a)') line(1:length)
              endif
           endif
        enddo
        write(*,"(a,23('*'),a,29('*'))") &
             '***', ' End of input data file '
        call io_close(lun_tmp)
!
!          "filein" for fdf is now the temporary file. 
!          This was necessary historically to allow
!          the rewinds involved in fdf operation.
!
     endif
  endif
! ...

! Set up fdf ...
!
! Choose a 'unique' prefix for the log (and possible debug) fdf files
! The 5-digit sequence might be slightly different in different
! processors, depending on the system time.
  call system_clock( count )
  write(fileout,"(a,i5.5,a)") 'fdf-', mod(count,100000), ".log"

  call fdf_init(filein,trim(fileout))

  ! Initialize the verbosity setting
  call init_verbosity('TBT.Verbosity',5)

! Define Name of the system ...
  sname = fdf_get('SystemName',' ')
  if (Node.eq.0) then
     write(*,'(/a,71("-"))') 'reinit: '
     write(*,'(a,a)') 'reinit: System Name: ',trim(sname)
     write(*,'(a,71("-"))') 'reinit: '
  endif
! ...

! Define System Label (short name to label files) ...
  slabel = fdf_get('SystemLabel','siesta')
  if (Node.eq.0) then
     write(*,'(a,a)') 'reinit: System Label: ',trim(slabel)
     write(*,'(a,71("-"))') 'reinit: '
  endif
! ...

end subroutine tbt_reinit
