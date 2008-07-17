! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine reinit(sname) 

C Subroutine to initialise the reading of the data for SIESTA 
C
C     It uses the FDF (Flexible Data Fromat) package 
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
      use precision
      use parallel,    only : Node
      use fdf
      use files,       only : slabel

      implicit none

      character(len=*) :: sname
      character string*20

C  Internal variables .................................................
      character
     .  filein*20, fileout*20

      character
     .  slabel_default*59, sname_default*20

      logical debug_input, file_exists

C Print Welcome and Presentation .......................................
      if (Node.eq.0) then
        write(6,'(/a)') 
     .   '                           ***********************       '
        write(6,'(a)') 
     .   '                           *  WELCOME TO SIESTA  *       '
        write(6,'(a)')
     .   '                           ***********************       '
C ..................
      endif

c
c     Choose proper file for fdf processing
c     (INPUT_DEBUG if exists or standard input)
c
      inquire(file='INPUT_DEBUG',exist=debug_input)
      if (debug_input) then
        if (Node.eq.0) then
          write(6,'(a)') 'WARNING: ' //
     $        'Siesta is reading its input from file INPUT_DEBUG'
        endif
        filein = 'INPUT_DEBUG'
      else
c
c       Read from standard input
c
        if (Node.eq.0) then
          write(6,'(/a)') 'reinit: Reading from standard input'
        endif
        filein = 'INPUT.fdf'
      endif
c ...

C Set up fdf ...
      fileout = 'out.fdf'
      call fdf_init(filein,fileout)

C Dump data file to output file
      if (Node.eq.0) then
        write(6,'(a,23(1h*),a,28(1h*))')
     .      '***', ' Dump of input data file '

!!!!!!!!! NEED TO BE DONE IN FDF LIBRARY !!!!!!!!!!!!!
!        call fdf_dumpAll()

        write(6,'(a,23(1h*),a,29(1h*))')
     .      '***', ' End of input data file '
      endif
C ...

C Define Name of the system ...
      sname_default = ' '
      sname = fdf_string('SystemName',sname_default)
      if (Node.eq.0) then
        write(6,'(/a,71(1h-))') 'reinit: '
        write(6,'(a,a)') 'reinit: System Name: ',trim(sname)
        write(6,'(a,71(1h-))') 'reinit: '
      endif
C ...

C Define System Label (short name to label files) ...
      slabel_default  = 'siesta'
      slabel = fdf_string('SystemLabel',slabel_default)
      if (Node.eq.0) then
        write(6,'(a,a)') 'reinit: System Label: ',slabel
        write(6,'(a,71(1h-))') 'reinit: '
      endif
C ...

      return
      end
