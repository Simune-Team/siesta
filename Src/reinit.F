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
      use siesta_master, only: input_file  ! fdf data file

      implicit none

      character(len=*), intent(out) :: sname

C  Internal variables .................................................
      character(len=20) filein, fileout

      integer  ::  count, length, lun, lun_tmp, iostat
      character :: line*256

      logical debug_input, file_exists

C Print Welcome and Presentation .......................................
!
      filein = "fdf_input"  

      if (Node.eq.0) then
        write(6,'(/a)') 
     .   '                           ***********************       '
        write(6,'(a)') 
     .   '                           *  WELCOME TO SIESTA  *       '
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
     $          'Siesta is reading its input from file INPUT_DEBUG'
           filein = 'INPUT_DEBUG'

        else

           if (input_file/=' ') then   ! Read data from given file
             write(6,'(a)') 'reinit: Reading from '//trim(input_file)
             call io_assign(lun)
             open(lun,file=input_file,form='formatted',status='old')
             rewind(lun)
           else                        ! Read from standard input
             write(6,'(/a)') 'reinit: Reading from standard input'
             lun = 5
           end if ! (input_file/=' ')

           call io_assign(lun_tmp)
           do  ! make sure we get a new file
              call system_clock( count )
              write(filein,"(a,i5.5)") 'INPUT_TMP.', mod(count,100000)
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
           
           if ( lun /= 5 ) then
              ! Be sure to close the input master file
              ! Otherwise we may end up opening too many files
              call io_close(lun)
           end if
c
!          "filein" for fdf is now the temporary file. 
!          This was necessary historically to allow
!          the rewinds involved in fdf operation.
!
        endif
      endif
c ...

C Set up fdf ...
!
!     Choose a 'unique' prefix for the log (and possible debug) fdf files
!     The 5-digit sequence might be slightly different in different
!     processors, depending on the system time.

      call system_clock( count )
      write(fileout,"(a,i5.5,a)") 'fdf-', mod(count,100000), ".log"

      call fdf_init(filein,trim(fileout))

C Define Name of the system ...
      sname = fdf_string('SystemName',' ')
      if (Node.eq.0) then
        write(6,'(/a,71("-"))') 'reinit: '
        write(6,'(a,a)') 'reinit: System Name: ',trim(sname)
        write(6,'(a,71("-"))') 'reinit: '
      endif
C ...

C Define System Label (short name to label files) ...
      slabel = fdf_string('SystemLabel','siesta')
      slabel = trim(slabel)
      ! Check that there are no spaces in the SystemLabel
      count = len_trim(slabel)
      length = index(slabel, ' ')
      if ( length > 0 .and. length < count ) then
         call die('SystemLabel must *NOT* contain any spaces!')
      end if
      if (Node.eq.0) then
        write(6,'(a,a)') 'reinit: System Label: ',slabel
        write(6,'(a,71("-"))') 'reinit: '
      endif
C ...

      end
