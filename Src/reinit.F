      subroutine reinit(slabel, sname) 
C *********************************************************************
C Subroutine to initialise the reading of the data for SIESTA 
C
C     It uses the FDF (Flexible Data Fromat) package 
C     of J.M.Soler and A.Garcia
C
C Taken from redata. Writen by P.Ordejon, December'96
C **************************** OUTPUT *********************************
C character*20 slabel      : System Label (to name output files)
C character*59 sname       : System Name
C **********************************************************************

C
C  Modules
C
      use precision
      use fdf
#ifdef MPI
      use mpi
#endif

      implicit none

      character
     .  slabel*20, sname*59

C  Internal variables .................................................
      character
     .  filein*20, fileout*20, line*150 

      character
     .  slabel_default*59, sname_default*20

      integer 
     .  Node, length, lun

      logical file_exists
#ifdef MPI
      integer
     .  MPIerror
#endif


C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

C Print Welcome and Presentation .......................................
      if (Node.eq.0) then
        write(6,'(/a)') 
     .   '                           ***********************       '
        write(6,'(a)') 
     .   '                           *  WELCOME TO SIESTA  *       '
        write(6,'(a)')
     .   '                           ***********************       '
C ..................

C Dump data file to output file .......................................
C and generate scratch file for FDF to read from
C

        call io_assign(lun)
        inquire(file='INPUT',exist=file_exists)
        if (.not. file_exists) then
           write(6,'(a)') 'WARNING: ' //
     $          'Siesta now reads its input from file INPUT...'
           call die('Input file not found')
        endif
        open(lun,file='INPUT',form='formatted',status='old')
        rewind(lun)

        write(6,'(/,a,18(1h*),a,28(1h*))')
     .    '*: ', ' Dump of input data file '
 10     continue
        read(lun,err=20,end=20,fmt='(a)') line
        call chrlen(line,0,length)
        if (length .ne. 0) write(6,'(a)') line(1:length)
        goto 10
 20     continue
        call io_close(lun)

        write(6,'(a,18(1h*),a,29(1h*))')
     .    '*: ', ' End of input data file '
C ..................

C Read data from FDF file..............................................

C Set up fdf ...
        filein = 'INPUT'
        fileout = 'out.fdf'
        call fdf_init(filein,fileout)
C ...

        write(6,'(/,a,18(1h*),a,30(1h*))')
     .    'redata: ', ' Simulation parameters '
        write(6,'(a)')  'redata:'
        write(6,'(a)')  'redata: The following are some of the '//
     .                           'parameters of the simulation.'
        write(6,'(a)')  'redata: A complete list of the parameters '//
     .                           'used, including defect values,'
        write(6,'(a,a)')'redata: can be found in file ',fileout
        write(6,'(a)')  'redata:'

C Defile Name of the system ...
        sname_default = ' '
        sname = fdf_string('SystemName',sname_default)
        write(6,'(a,71(1h-))') 'redata: '
        write(6,'(a,a)') 
     .   'redata: System Name: ',sname
        write(6,'(a,71(1h-))') 'redata: '
C ...

C Defile System Label (short name to label files) ...
        slabel_default  = 'siesta'
        slabel = fdf_string('SystemLabel',slabel_default)
        write(6,'(a,4x,a)') 
     .   'redata: System Label                     = ',slabel
C ...
      endif

C Global broadcast of strings
#ifdef MPI
      call MPI_Bcast(sname,59,MPI_character,0,MPI_Comm_World,
     .  MPIerror)
      call MPI_Bcast(slabel,20,MPI_character,0,MPI_Comm_World,
     .  MPIerror)
#endif

      return
      end

