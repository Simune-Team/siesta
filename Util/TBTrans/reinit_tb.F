! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine reinit(sname) 

C Subroutine to initialise the reading of the data for SIESTA 
C
C     It uses the FDF (Flexible Data Format) package 
C     of J.M.Soler and A.Garcia
C
C Taken from redata. Writen by P.Ordejon, December'96
C **************************** OUTPUT *********************************
C character    slabel      : System Label (to name output files)
C character(len=*) sname       : System Name
C **********************************************************************

C
C  Modules
C
      use parallel,    only : Node
      use fdf
      use files,       only : slabel

      implicit none

      character(len=*), intent(out) :: sname

C  Internal variables .................................................
      character(len=20) filein, fileout, string

      integer  ::  count, length, lun, lun_tmp, iostat
      character :: slabel_default*59, sname_default*20, line*256

      logical debug_input, file_exists

C Print Welcome and Presentation .......................................

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
        write(6,'(/a)') 
     .   '                           ***********************       '
        write(6,'(a)') 
     .   '                           *  WELCOME TO TBTRANS *       '
        write(6,'(a)')
     .   '                           ***********************       '
C ..................
!
!       Set name of file to read from. Done only
!       in the master node.
!
c
c     Choose proper file for fdf processing
c     (INPUT_DEBUG if it exists or "standard input",
c      processed and dumped to a temporary file)
c
        inquire(file='INPUT_DEBUG',exist=debug_input)
        if (debug_input) then
           write(6,'(a)') 'WARNING: ' //
     $          'TBTrans is reading its input from file INPUT_DEBUG'
           filein = 'INPUT_DEBUG'

        else
c
c          Read from standard input (dumped to a temp file)
c
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
           open(lun_tmp,file=filein,
     $          form='formatted',status='replace')
           rewind(lun_tmp)
           write(6,"(a,23('*'),a,28('*'))")
     .          '***', ' Dump of input data file '
c
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
           write(6,"(a,23('*'),a,29('*'))")
     .          '***', ' End of input data file '
           call io_close(lun_tmp)
c
!          "filein" for fdf is now the temporary file. 
!          This was necessary historically to allow
!          the rewinds involved in fdf operation.
!
        endif
      endif
c ...

C Set up fdf ...
      fileout = 'fdf.log'
      call fdf_init(filein,fileout)

C Define Name of the system ...
      sname_default = ' '
      sname = fdf_string('SystemName',sname_default)
      if (Node.eq.0) then
        write(6,'(/a,71("-"))') 'reinit: '
        write(6,'(a,a)') 'reinit: System Name: ',trim(sname)
        write(6,'(a,71("-"))') 'reinit: '
      endif
C ...

C Define System Label (short name to label files) ...
      slabel_default  = 'siesta'
      slabel = fdf_string('SystemLabel',slabel_default)
      if (Node.eq.0) then
        write(6,'(a,a)') 'reinit: System Label: ',slabel
        write(6,'(a,71("-"))') 'reinit: '
      endif

      end
