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
      use sys, only : bye

      implicit none

      character(len=*), intent(out) :: sname

C  Internal variables .................................................
      character(len=50) filein, fileout

      integer  ::  count, in, length, lun, lun_tmp, iostat
      character :: slabel_default*59, sname_default*20
      character(len=256) :: line, val


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

! ..................
!
!       Set name of file to read from. Done only
!       in the master node.
!
#ifndef NO_F2003
      count = command_argument_count()
      if ( count > 0 ) then
         filein = ' '
         call get_command_argument(count,filein,length)
         inquire(file=filein,exist=debug_input)
         if ( .not. debug_input ) count = 0
      end if
#endif

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

#ifndef NO_F2003
        else if ( count > 0 ) then

           ! Get file-name from input line
           filein = ' '
           call get_command_argument(count,filein,length)
           if ( length > len(filein) ) then
              call die ('The argument is too long to be retrieved, &
     &please limit your-self to 50 characters for the input file') 
           end if
           inquire(file=filein,exist=debug_input)
           if ( .not. debug_input ) then
              call die('Input file '//trim(filein)//' does not exist? &
     &Have you specified the wrong file-name?') 
           end if
           
           write(*,'(/,2a)') 'reinit: Reading from ',trim(filein)

#endif
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

#ifndef NO_F2003
! Read special variables from the command line
      count = command_argument_count()
      if ( count > 0 ) then
         in = 0
         do while ( in <= count - 1 )
            
            in = in + 1
            call get_command_argument(in,line,length)
            
            ! If it is not an option, skip it
            if ( line(1:1) /= '-' ) cycle

            do while ( line(1:1) == '-' )
               line = line(2:)
            end do

            ! We allow these line
            if ( line(1:3) == 'fdf' ) then
               if ( in >= count - 1 ) 
     &              call die('Missing argument on command line')
               in = in + 1
               call get_command_argument(in,line,length)
               
               ! We allow these variations:
               !  FDFLabel=0.1:eV
               !  FDFLabel:0.1:eV
               !  FDFLabel=0.1=eV
               line = cmd_tokenize(line)
               call fdf_overwrite(line)
#ifdef TRANSIESTA
            else if ( line(1:1) == 'V' ) then
               if ( in >= count - 1 )
     &              call die('Missing argument on command line')
               in = in + 1
               call get_command_argument(in,line,length)
               line = cmd_tokenize(line)
               line = 'TS.Voltage '//trim(line)
               call fdf_overwrite(line)
#endif
               
            else if ( line(1:1) == 'L' ) then
               if ( in >= count - 1 ) 
     &              call die('Missing argument on command line')
               in = in + 1
               call get_command_argument(in,line,length)
               line = cmd_tokenize(line)
               line = 'SystemLabel '//trim(line)
               call fdf_overwrite(line)
               
            else if ( line(1:4) == 'help' .or. line(1:1) == 'h' ) then
      write(*,'(a)')'Help for calling the tight-binding transport code'
      write(*,'(a)')'  -fdf <label>=<value>[:<unit>]'
      write(*,'(a)')'      Set the label to the corresponding value.'
#ifdef TRANSIESTA
      write(*,'(a)')'  -V <value>:<unit>'
      write(*,'(a)')'      Short-hand for setting TS.Voltage.'
#endif
      write(*,'(a)')'  -L <name>'
      write(*,'(a)')'      Short-hand for setting SystemLabel.'
      write(*,'(a)')'  <fdf-file>'
      write(*,'(a)')
     &     '      Use file as fdf-input, you need not to pipe it in.' 
      call bye('Help-menu requested, stopping')
            end if
            
         end do
      end if
#endif


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
C ...

      contains
      
      function cmd_tokenize(line) result(tline)
      character(len=*), intent(in) :: line
      character(len=len(line)) :: tline
      
      integer :: i, n
      n = len(tline)
      tline = line
      do i = 1 , n
         if ( tline(i:i) == ':' .or. 
     &        tline(i:i) == '=' ) then
            tline(i:i) = ' '
         end if
      end do
      end function cmd_tokenize
      end
