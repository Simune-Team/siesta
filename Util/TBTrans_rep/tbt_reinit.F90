! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! Transferred to be used in the TBTrans utility.
!
subroutine tbt_reinit(sname) 

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
  use files,       only : slabel

  implicit none

  character(len=*), intent(out) :: sname

!  Internal variables .................................................
  character(len=30) filein, fileout, string

  integer  ::  count, length, lun, lun_tmp, iostat
  character :: slabel_default*59, sname_default*20, line*256

  logical debug_input, file_exists

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
     write(6,'(a)') &
          '                           *  WELCOME TO TBTrans  * '
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
        write(6,'(a)') 'WARNING: ' // &
             'TBTrans is reading its input from file INPUT_DEBUG'
        filein = 'INPUT_DEBUG'

     else
!
!          Read from standard input (dumped to a temp file)
!
        write(6,'(/a)') 'reinit: Reading from standard input'
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
        write(6,"(a,23('*'),a,28('*'))") &
             '***', ' Dump of input data file '
!
        do
           read(lun,iostat=iostat,fmt='(a)') line
           if (iostat /= 0 ) exit
           length = len_trim(line)
           if (length /= 0) then
              write(6,'(a)') line(1:length)
              if (.not. debug_input) then
                 write(lun_tmp,'(a)') line(1:length)
              endif
           endif
        enddo
        write(6,"(a,23('*'),a,29('*'))") &
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
  write(fileout,"(a,i5,a)") 'fdf-', mod(count,100000), ".log"

  call fdf_init(filein,trim(fileout))

! Define Name of the system ...
  sname_default = ' '
  sname = fdf_string('SystemName',sname_default)
  if (Node.eq.0) then
     write(6,'(/a,71("-"))') 'reinit: '
     write(6,'(a,a)') 'reinit: System Name: ',trim(sname)
     write(6,'(a,71("-"))') 'reinit: '
  endif
! ...

! Define System Label (short name to label files) ...
  slabel_default  = 'siesta'
  slabel = fdf_string('SystemLabel',slabel_default)
  if (Node.eq.0) then
     write(6,'(a,a)') 'reinit: System Label: ',slabel
     write(6,'(a,71("-"))') 'reinit: '
  endif
! ...

end subroutine tbt_reinit
