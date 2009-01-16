      PROGRAM gengf

C ##################################################################
C ##                                                              ##       
C ##                 "Generate Greens Functions"                  ##
C ##                                                              ##       
C ##                                                              ## 
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##               Mikroelektronik Centret (MIC)                  ##
C ##           Technical University of Denmark (DTU)              ##
C ##                                                              ##
C ##################################################################
C
C Driver routine for generation of Greens functions on a complex
C contour.
C
C Copyright by Mads Brandbyge, 1999, 2000.
C The use of this program is allowed for not-for-profit research only.
C Copy or disemination of all or part of this package is not
C permitted without prior and explicit authorization by the author.
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C UNITS:
C     Ry, and Bohr is used throughout the program unless explicitly 
C     stated.
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use cont

c======================================================================
      implicit none
c======================================================================

      logical PRINTALOT
c      parameter(PRINTALOT=.FALSE.)
      parameter(PRINTALOT=.TRUE.)
      real*8 eV
      parameter ( eV = 1.d0/13.60580d0 )
      real*8 Pi
      parameter(Pi=3.141592653589793238462643383279503d0)


c======================================================================
c FROM FDF
      real*8 Volt               ! the bias voltage
      logical tcmplxcon, tbtjob

C OUT-FILE:
      integer joutfile          ! unit-number of out-file, OUT

c
c     Energy Contour:
c
      complex*16, dimension (:), pointer:: contour,wgf
c      complex*16 contour(NCONTOUR),WGF(NCONTOUR) ! Total "contour"

c     Bulk density of states:
      complex*16, dimension(:), allocatable :: zbulkdosL,zbulkdosR
c      complex*16 zbulkdosL(NCONTOUR), zbulkdosR(NCONTOUR)
c======================================================================
c Fermi levels
      real*8 EFermiL,EFermiR,EFermi0

c No. contour points (NCONTOUR if VOLT=0 else NCONTOUR + NVOLT)
      integer NEn

c ==================================================================
c Helpers, workspace, tempos, auxilaries etc...

      integer iEn
      integer Ncon
      logical isvolt

      integer ncontour,ncontour0,nvolt

      real*8 rtmp,totL,totR,totLonR,totRonL,totLra,totRra

c Files
      external io_assign, io_close
      integer unit2


c FDF-stuff:
      character*20 slabel       ! System Label (to name output files)
      character*59 sname        ! System Name
        
      integer length, lun
      character line*150, filein*20, fileout*20
      character  slabel_default*59, sname_default*20

      character paste*33,fname*33
      external paste
C Define FDF calls ....................................................
      include 'fdf/fdfdefs.h'

c=======================================================================
c BEGIN:
c=======================================================================

C Print Welcome and Presentation .......................................
      write(6,'(/a)') 
     . '                           **********************       '
      write(6,'(a)') 
     . '                           *  WELCOME TO GENGF  *       '
      write(6,'(a)')
     . '                           **********************       '
C ..................

c      call TSprversion

C Dump data file to output file .......................................
C and generate scratch file for FDF to read from
C
      write(6,'(/,a,18(1h*),a,28(1h*))')
     .  'gengf: ', ' Dump of input data file '

      call io_assign(lun)
      open(lun,file='FDF_STDIN',form='formatted',status='unknown')
      rewind(lun)

 10   continue
         read(5,err=20,end=20,fmt='(a)') line
         call chrlen(line,0,length)
         write(lun,'(a)') line(1:length)
         if (length .ne. 0) write(6,'(a)') line(1:length)
         goto 10
 20   continue
      call io_close(lun)

      write(6,'(a,18(1h*),a,29(1h*))')
     .  'gengf: ', ' End of input data file '
C ..................

C Read data from FDF file..............................................

C Set up fdf ...
      filein = 'FDF_STDIN'
      fileout = 'out.fdf'
      call fdf_init(filein,fileout)
C ...

C Defile Name of the system ...
      sname_default = 'genGF'
      sname = fdf_string('SystemName',sname_default)
      write(6,'(a,71(1h-))') 'gengf: '
      write(6,'(a,a)') 
     . 'gengf: System Name: ',sname
      write(6,'(a,71(1h-))') 'gengf: '
C ...

C Defile System Label (short name to label files) ...
      slabel_default  = 'gengf'
      slabel = fdf_string('SystemLabel',slabel_default)
      write(6,'(a,4x,a)') 
     . 'gengf: System Label                     = ',slabel
C ...


      joutfile=6                !write info to stdout
      
      Volt=fdf_double('TS.Voltage',0d0)
      write(6,'(a,f10.4,a)') 
     .     'gengf: Voltage  (Volt)           = ',Volt


c
c     Setup Energy contours       
c        
c     We take the common Fermi energy of the electrodes to be zero!
c     An initial equilibration between the electrodes is assumed.
c
         EFermiL =  Volt*eV/2.0d0 !Ry
         EFermiR = -Volt*eV/2.0d0 !Ry
         EFermi0 =       0.0d0
         isvolt=.false.
         if(dabs(Volt*eV).GT.0.001d0) isvolt=.true.

c
c     Setup contour
c
         call setupcontour(joutfile,NEn,
     .     EFermi0,EFermiL,EFermiR,
     .     WGF,contour,ncontour,ncontour0,nvolt)

         allocate(zbulkdosL(ncontour))
         allocate(zbulkdosR(ncontour))

      tcmplxcon = fdf_boolean('TS.ComplexContour',.true.)

      tbtjob = .not.tcmplxcon
c
c     Generate/Read-in Green's functions on contour
c
c -------------------------------------------------------
c Left GAAL:

      CALL green(joutfile,
     &     NEn,contour,WGF,EFermiL,zbulkdosL,.true., tbtjob )

c -------------------------------------------------------
c Right GAAR:

      CALL green(joutfile,
     &     NEn,contour,WGF,EFermiR,zbulkdosR,.false.,tbtjob )

c -------------------------------------------------------
c
c Integrate bulk dos on contours: (L,R)
c

      if(tcmplxcon) then
         ncon = NCONTOUR0 
      else
c         ncon = NCONTOUR
         ncon = nen
      endif

      if(tcmplxcon) then

c Left:          
      call io_assign(unit2)
      fname = paste(slabel,'.GFDOSLeft')
      OPEN(unit2,FILE=fname, STATUS='UNKNOWN')

      rtmp=0.0d0
      do iEn=1,ncon
         rtmp = rtmp - 
     &        (2d0/Pi)*DIMAG(zbulkdosL(iEn)*WGF(iEn))
         write(unit2,'(i3,f10.5,1X,f10.5,1X,f12.5,1X,f12.5)')
     &        iEn,rtmp,contour(iEn),(2d0/Pi)*DIMAG(zbulkdosL(iEn))

      end do
      totL=rtmp
      call io_close(unit2)

c Right:
      call io_assign(unit2)
      fname = paste(slabel,'.GFDOSRight')
      OPEN(unit2,FILE=fname, STATUS='UNKNOWN')

      rtmp=0.0d0
      if(isvolt) then
         do iEn=NCONTOUR0+1,2*NCONTOUR0
            rtmp = rtmp - 
     &        (2d0/Pi)*DIMAG(zbulkdosR(iEn)*WGF(iEn))
            write(unit2,'(i3,f10.5,1X,f10.5,1X,f12.5,1X,f12.5)')
     &           iEn,rtmp,contour(iEn),(2d0/Pi)*DIMAG(zbulkdosR(iEn))
         end do
      else
         do iEn=1,ncon
            rtmp = rtmp - 
     &        (2d0/Pi)*DIMAG(zbulkdosR(iEn)*WGF(iEn))
            write(unit2,'(i3,f10.5,1X,f10.5,1X,f12.5,1X,f12.5)')
     &           iEn,rtmp,contour(iEn),(2d0/Pi)*DIMAG(zbulkdosR(iEn))
         end do
      end if                    !volt>0
      totR=rtmp
      call io_close(unit2)

      write(joutfile,*) 'Total charge: L,R:', totL,totR


c Non-equilibrium check:
      if(isvolt) then
         
         call io_assign(unit2)
         fname = paste(slabel,'.GFDOSNeq')
         OPEN(unit2,FILE=fname, STATUS='UNKNOWN')  
         
c     L-lead on R contour:
         rtmp=0.0d0
         do iEn=NCONTOUR0+1,2*NCONTOUR0
            rtmp = rtmp - 
     &           (2d0/Pi)*DIMAG(zbulkdosL(iEn)*WGF(iEn))
         end do
         totLonR=rtmp

c     L-lead on real-axis from EFR to EFL:
         rtmp=0.0d0
         write(unit2,*) 'Left real-axis:'
         do iEn=2*NCONTOUR0+1,2*NCONTOUR0+NVOLT
            rtmp = rtmp - 
     &           (2d0/Pi)*DIMAG(zbulkdosL(iEn)*WGF(iEn))
            write(unit2,'(i3,f10.5,1X,f10.5,1X,f12.5,1X,f12.5)')
     &     iEn,rtmp,contour(iEn),(2d0/Pi)*DIMAG(zbulkdosL(iEn))
         end do
         totLra=rtmp

c     R-lead on L contour:
         rtmp=0.0d0
         do iEn=1,NCONTOUR0
            rtmp = rtmp - 
     &           (2d0/Pi)*DIMAG(zbulkdosR(iEn)*WGF(iEn))
         end do
         totRonL=rtmp

c     R-lead on real-axis from EFR to EFL:
         write(unit2,*) 'Right real-axis:'
         rtmp=0.0d0
         do iEn=2*NCONTOUR0+1,2*NCONTOUR0+NVOLT
            rtmp = rtmp - 
     &           (2d0/Pi)*DIMAG(zbulkdosR(iEn)*WGF(iEn))
            write(unit2,'(i3,f10.5,1X,f10.5,1X,f12.5,1X,f12.5)')
     &      iEn,rtmp,contour(iEn),(2d0/Pi)*DIMAG(zbulkdosR(iEn))
         end do
         totRra=rtmp

         write(joutfile,*) ' '         
         write(joutfile,*) 'Charge in voltage window: '
         write(joutfile,*) 'ContourL-ContourR, Real-axis integr.:'
         write(joutfile,*) 'LEFT: ', totL-totLonR,totLra
         write(joutfile,*) 'RIGHT:', totRonL-totR,totRra
         write(joutfile,*) 'Estimated real-axis integration error:'
         rtmp = 100d0*dabs(1d0-(totLra/(totL-totLonR)))
         write(joutfile,'(a7,F6.2,a2)') 'LEFT:  ',rtmp,' %' 
         rtmp = 100d0*dabs(1d0-(totRra/(totRonL-totR)))
         write(joutfile,'(a7,F6.2,a2)') 'RIGHT: ',rtmp,' %' 

         call io_close(unit2) 
        
      end if                    !volt

       else ! tcmplxcon



c Left:
      call io_assign(unit2)
      fname = paste(slabel,'.GFDOSLeft')
      OPEN(unit2,FILE=fname, STATUS='UNKNOWN')

      rtmp=0.0d0
      do iEn=1,ncon
         rtmp = rtmp -
     &        (2d0/Pi)*DIMAG(zbulkdosL(iEn)*WGF(iEn))
         write(unit2,'(i3,f10.5,1X,f10.5,1X,f12.5,1X,f12.5)')
     &        iEn,rtmp,contour(iEn),(2d0/Pi)*DIMAG(zbulkdosL(iEn))

      end do
      totL=rtmp
      call io_close(unit2)

c Right:
      call io_assign(unit2)
      fname = paste(slabel,'.GFDOSRight')
      OPEN(unit2,FILE=fname, STATUS='UNKNOWN')

      rtmp=0.0d0
         do iEn=1,ncon
            rtmp = rtmp -
     &        (2d0/Pi)*DIMAG(zbulkdosR(iEn)*WGF(iEn))
            write(unit2,'(i3,f10.5,1X,f10.5,1X,f12.5,1X,f12.5)')
     &           iEn,rtmp,contour(iEn),(2d0/Pi)*DIMAG(zbulkdosR(iEn))
         end do
      totR=rtmp
      call io_close(unit2)

      write(joutfile,*) 'Total charge: L,R:', totL,totR


      endif










c     
c -------------------------------------------------------
c

c      call io_close(joutfile)

c ===================================================================
      end
c ===================================================================






