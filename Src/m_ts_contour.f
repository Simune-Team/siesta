      MODULE m_ts_contour
!
! Routines that are used to setup the contour for integration of the GFs
!
!=============================================================================
! CONTAINS:
!          1) distrcontour
!          2) mkCplxContour
!          3) mkRealContour
!          4) setupcontour


      implicit none

      public :: setupcontour, distrcontour, mkRealContour

      private

      CONTAINS



C ##################################################################
C ##                                                              ## 
C ##   Distribute contour points on parallel nodes                ##
C ##                                                              ##   
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ##################################################################

      subroutine distrcontour(
     .     WGF,contour,NCONTOUR,NCONTOUR0,NVOLT,isvolt,Nodes,
     .     NPARACONTOUR,paracontour,paraWGF,
     .     contourpart,points2read)

       use parallel, only : Node

c======================================================================
      use precision, only : dp
      implicit none
c======================================================================
c     The main idea is to "fill" up points with no weight on excess processors. 
c     "Easy" equilibrium points are done first on all processors and 
c     the "Hard" non-equilibrium points are done secondly on all processors.
c======================================================================      

      logical PRINTALOT
      parameter (PRINTALOT=.True.)


c     INPUT

c     Simple Energy Contour and weight:

      integer NCONTOUR          ! No. of "serial" contour points
      integer NVOLT             ! No. voltage points
      integer NCONTOUR0         ! No. of points on equilibrium contour
      logical isvolt            ! True if finite voltage
      integer Nodes             ! No. nodes for parallel processing
      complex(dp) WGF(NCONTOUR), contour(NCONTOUR)


c======================================================================
c     OUTPUT
      integer    NPARACONTOUR   ! No. points on parallel contour
      complex(dp) paracontour(1:NCONTOUR,0:Nodes-1) ! Energy points distributed on nodes
      complex(dp) paraWGF(1:NCONTOUR,0:Nodes-1)     ! Weights on distributed energypoints 
      character contourpart(1:NCONTOUR,0:Nodes-1)  ! 'L' if belonging to Left Eq. contour-part
c                                                  ! 'R' if belonging to Right Eq. contour-part
c                                                  ! 'V' if belonging to Non-eq. contour-part      
c                                                  ! 'N' if no voltage present
      integer points2read(1:NCONTOUR)              ! No. points to read-in from GF-file for
c                                                    the actual parallel contour point.
c======================================================================
c     Helpers
      integer i,ipe,iEn,itot,inode

      integer itotp2r,itotr
c=======================================================================
c BEGIN:
c=======================================================================

      call timer('dcontour',1)

c Init

      points2read = 0
      contourpart = 'X'


c ------------------ no voltage ------------------------
      if(.not. isvolt) then
         
         iEn=0                  !actual contour points (in GF file)
         ipe=0                  !energy points in parallel
         itot=0                 !total no. points done in all

         do i=1,ceiling(1.d0*NCONTOUR0/Nodes)
            ipe = ipe+1 
           
            do inode=0,Nodes-1
               
               itot=itot+1      
               if(itot .le. NCONTOUR0) then
                  iEn = iEn+1
                  paraWGF(ipe,inode) = WGF(iEn)
                  paracontour(ipe,inode) = contour(iEn)
                  points2read(ipe) = points2read(ipe) + 1
               else
                  paraWGF(ipe,inode) = dcmplx(0d0,0d0) !points with no weight!!
                  paracontour(ipe,inode) = contour(NCONTOUR0)
               end if
               
               contourpart(ipe,inode) = 'N'
      
            end do              ! Nodes
         end do                 ! parallel contourpoints
       

      end if                    ! no voltage
c ------------------ no voltage END ---------------------


      if(Node.eq.0) write(6,*) NCONTOUR,NCONTOUR0

c ------------------ Finite Voltage ---------------------
         
      if(isvolt) then
         
         iEn=0                  !actual contour points (in GF file)
         ipe=0                  !energy points in parallel
         itot=0                 !total no. points done in all
         

c     Distribute the "easy" Left/Right points on nodes

         do i=1,ceiling(2d0*NCONTOUR0/Nodes) 
            ipe=ipe+1           
            
            do inode=0,Nodes-1
               
               itot=itot+1                     
               if(itot .le. 2*NCONTOUR0) then
                  iEn = iEn+1
                  paraWGF(ipe,inode) = WGF(iEn)
                  paracontour(ipe,inode) = contour(iEn)
                  points2read(ipe)=points2read(ipe)+1
               else
                  paraWGF(ipe,inode) = dcmplx(0d0,0d0) !points with no weight!!
                  paracontour(ipe,inode) = contour(2*NCONTOUR0)
               end if

               contourpart(ipe,inode) = 'R'
               if(itot .le. NCONTOUR0)  contourpart(ipe,inode) = 'L'

            end do              ! Nodes
         end do                 ! parallel contourpoints
         


c     Distribute the "hard" points on nodes 


            itotr=iEn     ! no. actual contour points so far

         do i=1,ceiling(1.d0*NVOLT/Nodes) 
            ipe=ipe+1           
            
            do inode=0,Nodes-1
               itot=itot+1      
               itotr=itotr+1      
               if(itotr .le. 2*NCONTOUR0 + NVOLT) then
                  iEn = iEn+1
                  paraWGF(ipe,inode) = WGF(iEn)
                  paracontour(ipe,inode) = contour(iEn)
                  points2read(ipe)=points2read(ipe)+1
               else
                  paraWGF(ipe,inode) = dcmplx(0d0,0d0) !points with no weight!!
                  paracontour(ipe,inode) = contour(2*NCONTOUR0 + NVOLT)
               end if
               
               contourpart(ipe,inode) = 'V'

            end do              ! Nodes
         end do                 ! parallel contourpoints
         
      end if                    ! voltage!

c ------------------ Finite Voltage END ------------------



      NPARACONTOUR = ipe        !No. parallel contour points

      if(PRINTALOT) then
 
         if(Node.eq.0) then
         write(*,*) 
     .  '----- DISTRIBUTION OF ENERGY POINTS AMONG PROCESSORS -----'
         write(*,*) 
     .        ' Node',' Point ',' Part',' Read-ins ',
     .        '        ZEnergy   ','              Weight' 


         itotp2r =0

         do ipe=1,NPARACONTOUR              
                 
               itotp2r = itotp2r + points2read(ipe)
               do inode=0,Nodes-1

               write(*,
     . '(I5,1X,I5,5X,A,2X,i5,5X,F9.5,1X,F9.5,5X,F9.5,1X,F9.5)')
     .              inode,ipe,contourpart(ipe,inode),points2read(ipe),
     .              paracontour(ipe,inode), paraWGF(ipe,inode)
               
            end do
         end do
         write(*,*) 'Total no. points:  ',itot
         write(*,*) 'Total no. points to read:  ',itotp2r,iEn
         write(*,*) 
     . '----------------------------------------------------------'

      end if
      end if

      call timer('dcontour',2)

c ===================================================================
      return
      end subroutine distrcontour
c ===================================================================


!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------



C ##################################################################
C ##      Generate complex energy contour for integration of      ##
C ##        density-matrix from retarded Greens function          ##
C ## For finite Voltage the non-equilibrium parts are integrated  ##
C ## (close-to) real axis.                                        ##
C ##      Finite temp. treated with the Sommerfeld expansion      ##
C ##      or using gauss fermi contours in the ends               ##   
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##         and  Kurt Stokbro                                    ##     
C ##################################################################

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Equilibrium:
c Modified Hans Skriver contour. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine mkCplxContour(itype,E1,E2,
     .     kT,GFeta,nencont,zcontour,wgzcontour,
     .     ncircle,nline,npol)

      use parallel, only : Node
      use sys, only : die
      use units, only : Pi
      use precision, only : dp
      use m_ts_aux_rout, only : nf, nf1, gaufermi0, gaufermi2, 
     .                           gaufermi10, gaufermi20, gauss

      implicit none
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Modified Hans Skriver:
      integer NT
      parameter(NT=10)          ! start line in modified HS at E2-NT*kT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Modified REAL AXIS CONTOUR:
      integer NGAUF,ngau,NTGAU
      parameter(NGAUF=8)          ! number of points [-inf,E2+NT*kT]
      parameter(NTGAU=2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer itype             !Contour type
c itype=1: Modified Hans Skriver 
c itype=2: 1. order Sommerfeld expansion 
c itype=3: gauss fermi integral in the ends
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INPUT        
      integer nencont           ! No. contour points
      integer ncircle,nline,npol
      real(dp) E1,E2              ! energy parameters 
      real(dp) kT                 ! temperature in Ry
      real(dp) GFeta              ! state broadening in Ry

c OUTPUT
      complex(dp) zcontour(nencont) ! points for GF
      complex(dp) wgzcontour(nencont) ! weights on GF
      
c     Helpers
      complex(dp), dimension(:), allocatable :: zc,wc,zl,wl,zp,wp
c      complex*16 zc(NCIRCLE), wc(NCIRCLE)
c      complex*16 zl(NLINE), wl(NLINE) 
c      complex*16 zp(NPOL), wp(NPOL)
      
      real(dp) delta,gamma,D
c      real*8 theta(NCIRCLE),x(NLINE)
c      real*8 wt(NCIRCLE+NLINE)
      real(dp) R, alpha,beta
      real(dp) wlt(NGAUF),xlt(NGAUF)

      
      real(dp), dimension(:), allocatable :: theta,x,wt
      complex(dp) ztmp, z0
      real(dp) min,max,rtmp,etaSF
      real(dp) EE1,EE2

      integer job,i,ic,j,Ni,ntgauuse


c=================================================================
c     BEGIN
c=================================================================
      call timer('CXcontour',1)
C Get Node number

! TODO... there is a lot of writing to temporary arrays
! This is not needed... NPA




 666  format(a10,f12.5,1X,f12.5,2X,f15.9,1X,f15.9,i4)

c------------------------------------------------------------
c
c     Modified Hans Skriver Contour
c
      if(itype.EQ.1) then       ! modified HS

         if(nencont.ne.NCIRCLE+NLINE+NPOL) then
            if(Node.eq.0)
     .      write(*,*) 
     &           'ERROR: MKCONTOUR:  nencont=',nencont,
     &           ', NCIRCLE+NLINE+NPOL=',NCIRCLE+NLINE+NPOL
            call die( 'ERROR: MKCONTOUR nencont not OK!' )
         end if

c               
c     Parameters
c     
         D = E2-E1
         Delta=NPOL*2.0d0*Pi*kT
         gamma=NT*kT

         alpha=dATAN(Delta/(D-gamma))
         
         R = dsqrt(Delta*Delta + (D - gamma)*(D - gamma))/ 
     &        (2d0*Cos(alpha))
         
         z0 = dcmplx(E1 + R, 0d0)
         beta=dasin(Delta/R)
           
         ic=0                   !contour index

c
c     Residuals:
c        
           if(Node.eq.0) 
     .   write(*,*) 'contour:  Residuals: '
         
         allocate(zp(npol))
!        call memory('A','Z',2*npol,'mkCplxContour')
         allocate(wp(npol))
         call memory('A','Z',2*npol,'mkCplxContour')

         do i=1,NPOL 
            zp(i)=dcmplx(E2,Pi*kT*(2.0d0*(i-1)+1d0))
            wp(i)=dcmplx(0d0,2d0*Pi*kT) 
            ic=ic+1
            zcontour(ic)=zp(i)
            wgzcontour(ic)=wp(i)           

           if(Node.eq.0) 
     .      write(*,666) 'contour: ',zp(i),wp(i),i

         end do                 !i
         call memory('D','Z',2*npol,'mkCplxContour')
         deallocate(zp)
         deallocate(wp)

c
c     Line contour:
c        
         if(Node.eq.0) 
     .   write(*,*) 'contour:  Fermi Line: '
         allocate (wt(nline))
         allocate (x(nline))
         call memory('A','D',2*nline,'mkCplxContour')         
         if(NT.EQ.10) then
            call gaufermi10(NLINE,x,wt)
         elseif(NT.EQ.20) then
            call gaufermi20(NLINE,x,wt)
         else
            if(Node.eq.0) then
            write(*,*) 'ERROR: ' 
            write(*,*) 
     &           'No Gauss quadrature for Fermi function '
            endif
            call die('No Gauss quadrature for Fermi function ')
         end if
         allocate(zl(nline))
         allocate(wl(nline))
         call memory('A','Z',2*nline,'mkCplxContour')         
         do i=1,NLINE
            j=NLINE-i+1         !reverse
            zl(i) = dcmplx(x(j)*kT + E2,Delta)
            wl(i) = -wt(j)*kT*dcmplx(1d0,0d0)
            ic=ic+1
            zcontour(ic)=zl(i)
            wgzcontour(ic)=wl(i)
            if(Node.eq.0)
     .      write(*,666) 'contour: ',zl(i),wl(i),i
         end do                 !ia
         call memory('D','D',2*nline,'mkCplxContour')         
         deallocate(wt) 
         deallocate(x)
         call memory('D','Z',2*nline,'mkCplxContour')         
         deallocate(zl)
         deallocate(wl)
         
c
c     Circle contour:
c
         if(Node.eq.0) 
     .   write(*,*) 'contour:  Circle: '
         min = beta 
         max = Pi
         job=0
         allocate (wt(ncircle))
         allocate (theta(ncircle))
         call memory('A','D',2*ncircle,'mkCplxContour') 
         call gauss(NCIRCLE, job, min, max, theta, wt)       

         allocate(zc(ncircle))
         allocate(wc(ncircle))
         call memory('A','Z',2*ncircle,'mkCplxContour') 
         do i=1,NCIRCLE 
            j=i
            ztmp=R*exp(theta(j)*dcmplx(0d0,1d0))
            zc(i) = z0 + ztmp
            wc(i)=wt(j)*nf((zc(i)-E2)/(kT))*dcmplx(0d0,1d0)*ztmp
            ic=ic+1
            zcontour(ic)=zc(i)
            wgzcontour(ic)=wc(i)

            if(Node.eq.0) 
     .      write(*,666) 'contour: ',zc(i),wc(i),i

         end do                 !i
         call memory('D','D',2*ncircle,'mkCplxContour') 
         call memory('D','Z',2*ncircle,'mkCplxContour') 
         deallocate(wt)
         deallocate(theta)
         deallocate(zc)
         deallocate(wc)
c------------------------------------------------------------
c
c     1. order Sommerfeld expansion using 2kT increments
c      
c We assume that the function to be integrated varies slowly on the
c kT-scale         

      else if(itype.EQ.2) then  ! 1. order Sommerfeld
         
         if(nencont.le.4) then
           if(Node.eq.0)
     .     write(6,*) 
     .     'ERROR: No. points=',nencont,' not valid for Sommerfeld'
           call die('ERROR: Contour: no. points not OK for Sommerfeld') 
         else
            Ni = nencont-2

            if(E1.GT.E2) then
               EE2=E1
               EE1=E2
            else
               EE2=E2
               EE1=E1
            end if

            delta = (EE2 - EE1)/(Ni-1) 
            
            if(GFeta .le. 0d0) call die(' ERROR: GFeta <= 0 ')

            etaSF = (kT)

            rtmp = (kT)*(kT)*Pi*Pi/(12.d0*etaSF)

            zcontour(1)     = dcmplx(EE1 - etaSF,GFeta)
            wgzcontour(1)   = dcmplx(0.25d0*delta + rtmp,0d0)

            zcontour(2)     = dcmplx(EE1 + etaSF,GFeta)
            wgzcontour(2)   = dcmplx(0.25d0*delta - rtmp,0d0)

            zcontour(Ni+1)  = dcmplx(EE2 - etaSF,GFeta)
            wgzcontour(Ni+1)= dcmplx(0.25d0*delta - rtmp,0d0)

            zcontour(Ni+2)  = dcmplx(EE2 + etaSF,GFeta)
            wgzcontour(Ni+2)= dcmplx(0.25d0*delta + rtmp,0d0)

            do i=3,Ni
               zcontour(i)   = dcmplx(delta*(i-2) + EE1,GFeta) 
               wgzcontour(i) = dcmplx(delta,0d0)             
            end do              !i
         end if                 ! nencont >=4

         if(E1.GT.E2) then
            do i=1,Ni+2
               wgzcontour(i)= -wgzcontour(i)
            end do
         end if                 !E1>E2


         if(Node.eq.0) then
         write(*,*) 'contour:  Non-equilibrium: '
         do i=1,Ni+2
            write(*,666) 'contour: ',zcontour(i),wgzcontour(i),i
         end do
         endif


c------------------------------------------------------------
c
c     Gaussian quadrature, using gaufermi in the ends
c      
c We assume that the function to be integrated varies slowly on the
c kT-scale         

      else if(itype.EQ.3) then  ! gaussian quadrature
c first determine how many points to use for the fermiline
         ntgauuse=NTGAU
         if (dabs(E2 - E1) .le. 3.d0*NTGAU *kT) then
            ntgauuse=0
         endif
         ngau=((2.d0+ntgauuse)* kT/dabs(E1-E2))*nencont +1.0d0
         if(ngau .gt. NGAUF+ntgauuse) ngau=NGAUF+ntgauuse
         if(nencont.lt.2*ngau+8) then
            ngau=(nencont-8)/2
         endif
         if (Node.eq.0)
     &   write(*,*) 'Gaussian fermi contour, ngau=',ngau
         if(ngau .le. 0) then
           if (Node.eq.0) write(*,*) 
     &       'ERROR: No. points=',nencont,
     &       ' not valid for real axis gaussian quadrature, min',
     &        10
          call die('ERROR: Contour: too few points for real axis int') 
         else
            Ni = nencont-ngau

            if(E1.GT.E2) then
               EE2=E1
               EE1=E2
            else
               EE2=E2
               EE1=E1
            end if

            if(GFeta .le. 0d0) call die(' ERROR: GFeta <= 0 ')

            if(ntgauuse .EQ.2) then
               call gaufermi2(ngau,xlt,wlt)
            elseif(ntgauuse .EQ. 0) then
               call gaufermi0(ngau,xlt,wlt)
            else
               if(Node.eq.0) then
               write(*,*) 'ERROR: ' 
               write(*,*) 
     &              'No Gauss quadrature for Fermi function '
               endif
               call die ('No Gauss quadrature for Fermi function')
            end if

            do i=1,ngau
               j=ngau-i+1      !reverse
               zcontour(i) = dcmplx(-xlt(j)*kT + EE1,GFeta)
               wgzcontour(i) = wlt(j)*kT*dcmplx(1d0,0d0)
               zcontour(Ni+ngau+1-i) = dcmplx(xlt(j)*kT + EE2,GFeta)
               wgzcontour(Ni+ngau+1-i) = wlt(j)*kT*dcmplx(1d0,0d0)
            enddo
            gamma=ntgauuse*kT
c set boundaries for gaussian quadrature

            delta = (EE2 - EE1-2.*gamma)/(Ni-ngau-1)
            do i=ngau+1,Ni
               rtmp=delta*(i-ngau-1) + EE1+gamma
               zcontour(i)   = dcmplx(rtmp,GFeta) 
               wgzcontour(i) = dcmplx(delta,0d0)*
     &              (nf1((rtmp-EE2)/kT) -nf1((rtmp-EE1)/kT))
            end do              !i
c extended simpsons rule
            wgzcontour(ngau+1)=wgzcontour(ngau+1)*17.d0/48.d0
            wgzcontour(ngau+2)=wgzcontour(ngau+2)*59.d0/48.d0
            wgzcontour(ngau+3)=wgzcontour(ngau+3)*43.d0/48.d0
            wgzcontour(ngau+4)=wgzcontour(ngau+4)*49.d0/48.d0
            wgzcontour(Ni+1-1)=wgzcontour(Ni+1-1)*17.d0/48.d0
            wgzcontour(Ni+1-2)=wgzcontour(Ni+1-2)*59.d0/48.d0
            wgzcontour(Ni+1-3)=wgzcontour(Ni+1-3)*43.d0/48.d0
            wgzcontour(Ni+1-4)=wgzcontour(Ni+1-4)*49.d0/48.d0
         end if                 ! nencont >=31

         if(E1.GT.E2) then
            do i=1,Ni+ngau
               wgzcontour(i)= -wgzcontour(i)
            end do
         end if                 !E1>E2

         if(Node.eq.0) then
         write(*,*) 'contour:  Non-equilibrium: '
         do i=1,Ni+ngau
            write(*,666) 'contour: ',zcontour(i),wgzcontour(i),i
         end do
         end if
c
c     itype not appropriate:
c
      else
         if(Node.eq.0)
     .   write(6,*) 'ERROR: mkCplxContour: Contour not defined'
         call die('ERROR:  mkCplxContour: Contour not defined') 
      end if
      

      call timer('CXcontour',2)        
C =========================================================        
      return
      end subroutine mkCplxContour
C =========================================================

!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


C ##################################################################
C ##      Generate (close to) real axis energy contour            ##
C ##           Transmission and DOS calculation                   ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################


      subroutine mkRealContour(E1,E2,
     &     GFeta,nencont,zcontour,wgzcontour)


      use parallel, only : Node
      use precision, only: dp
      use units, only : eV

      implicit none
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT 
      integer nencont           ! No. contour points
      real(dp) E1,E2              ! energy parameters 
      real(dp) GFeta              ! state broadening in Ry


c OUTPUT
      complex(dp) zcontour(nencont) ! points for GF
      complex(dp) wgzcontour(nencont) ! weights on GF
      
      real(dp) delta

      integer ic


c=======================================================================
c BEGIN:
c=======================================================================
      call timer('Rcontour',1)
C Get Node number


 666  format(a10,f12.5,1X,f12.5,2X,f15.9,1X,f15.9,i4)


c
c     Simple Line
c
    
      delta = (E2-E1)/(1.d0*max((nencont-1),1))
      write(*,*) 'contour:  Simple Line Contour:'
        
      do ic=1,nencont
         zcontour(ic) = dcmplx(E1+(ic-1)*delta,GFeta)
         wgzcontour(ic) = dcmplx(delta,0d0)
         if(Node.eq.0) then
         write(*,666) 'contour: ',
     &        dreal(zcontour(ic)),dimag(zcontour(ic)),
     &        dreal(wgzcontour(ic)),dimag(wgzcontour(ic)),ic
         
         endif
      end do                    !ic

      call timer('Rcontour',2)
C =========================================================
      return
      end subroutine mkRealContour
C =========================================================

!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------



C ##################################################################
C ##                                                              ## 
C ##   Setup complex energy contour by calling contour routines   ##
C ##                                                              ##   
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ## Changed to F90 by Jose-Luis Mozos, jlm@kanigo.icmab.es       ##
C ##################################################################

      subroutine setupcontour(NEn,
     &     EFermi0,EFermiL,EFermiR,
     &     WGF,contour,ncontour,ncontour0,
     &     nline,ncircle,npol,nvolt,
     &     Emin,GFEta,kT,smethod)

      use precision, only : dp
      use files, only : slabel
      use parallel, only : Node
      use units, only : eV

c======================================================================
      implicit none
c======================================================================
c INPUT
      real(dp) EFermi0            !Eq. Fermi energy
      real(dp) EFermiL,EFermiR    !Left/Right Fermi energies (finite voltage)
      
      integer ncontour, ncontour0
      integer nline,ncircle, npol,nvolt
      real(dp) Emin, GFEta, kT
      character(20) smethod
c OUTPUT

      complex(dp), dimension (:), pointer:: contour,wgf ! Total "contour"

      integer NEn               ! no. contour points to be used 
c                               ! (may differ from NCONTOUR)



c Helpers etc.
      logical isvolt            !True if finite voltage

      


      integer imethod           ! Complex contour method
      


      real(dp) Volt
C      real*8 temp_default       !default temperature
C      parameter ( temp_default  = 1.900d-3 ) 


c      complex*16 contmp(NCONTOUR),WGFtmp(NCONTOUR) 
      complex(dp), dimension (:), allocatable:: contmp,wgftmp 

      integer i,jc,ju
      
      external io_assign,io_close

      character paste*33,fname*33
C

c=======================================================================
c BEGIN:
c=======================================================================

      call timer('SUcontour',1)
C Get Node number


      isvolt=.false.
      Volt=EFermiL-EFermiR
      if(dabs(Volt).GT.0.001d0) isvolt=.true.


        if(smethod .eq. 'gaussfermi') then 
           imethod = 3
        else if(smethod .eq. 'sommerfeld') then
           imethod = 2 
        end if
  

c
c complex contour for integration of DOS
c

        if(isvolt) then
           ncontour=2*(npol+nline+ncircle)+nvolt
        else
           ncontour=  (npol+nline+ncircle)
        end if                 !volt
        ncontour0 = npol+nline+ncircle

        allocate(contour(ncontour)) 
        allocate(wgf(ncontour))
        call memory('A','Z',2*ncontour,'setupcontour')          

c
c Equilibrium Contour
c
        if(.not.isvolt) then
           allocate(contmp(ncontour0))
           allocate(wgftmp(ncontour0))
           call memory('A','Z',2*ncontour0,'setupcontour') 
! TODO Ask mads about this factor of eV...
           call mkCplxContour(1,Emin-Volt,EFermi0, ! removed factor of Volt*eV (EFermi[L,R] are in Ry)
     &           kT,GFeta,NCONTOUR0,contmp,WGFtmp, !Equilibrium part
     &           ncircle,nline,npol)
           jc=0
           do i=1,NCONTOUR0
              jc=jc+1
              contour(jc)=contmp(i)
c     
c     Note we put a minus here because the integral we want is the
c     negative of the pole-sum and C+L integral!!
c     
              WGF(jc)=-WGFtmp(i)
               
           end do              !i
           call memory('D','Z',2*ncontour0,'setupcontour') 
           deallocate(contmp)
           deallocate(wgftmp)

        else                   !finite voltage
c
c Non-Equilibrium Contour
c 

c Contour sequence:  
c (i)   equilibrium contour corresponding to LEFT EFermiL 
c (ii)  equilibrium contour corresponding to RIGHT EFermiR 
c (iii) real-axis contour between EFermiL and EfermiR

c (i):
           if(Node.eq.0) write(*,*) 'Left equilibrium:'

           allocate(contmp(ncontour0))
           allocate(wgftmp(ncontour0))
           call memory('A','Z',2*ncontour0,'setupcontour') 
           call mkCplxContour(1,
     &          Emin + (EFermiL-EFermi0),EFermiL,
     &          kT,GFeta,NCONTOUR0,contmp,WGFtmp,
     &          ncircle,nline,npol)
           jc=0
           do i=1,NCONTOUR0
              jc=jc+1
              contour(jc)=contmp(i)
c     
c     Note we put a minus here because the integral we want is the
c     negative of the pole-sum and C+L integral!!
c     
              WGF(jc)=-WGFtmp(i)
           end do              !i

c (ii):
           if(Node.eq.0) write(*,*) 'Right equilibrium:'

           call mkCplxContour(1,
     &          Emin + (EFermiR-EFermi0),EFermiR,
     &          kT,GFeta,NCONTOUR0,contmp,WGFtmp,
     &          ncircle,nline,npol)
           do i=1,NCONTOUR0
              jc=jc+1
              contour(jc)=contmp(i)
c     
c     Note we put a minus here because the integral we want is the
c     negative of the pole-sum and C+L integral!!
c     
              WGF(jc)=-WGFtmp(i)
           end do              !i
           call memory('D','Z',2*ncontour0,'setupcontour') 
           deallocate(contmp)
           deallocate(wgftmp)
           call memory('A','Z',2*nvolt,'setupcontour') 
           allocate(contmp(nvolt))
           allocate(wgftmp(nvolt))

c (iii):
           if(Node.eq.0) 
     .         write(*,*) 'Real-axis: From EFermiR to EFermiL'

         
           call mkCplxContour(imethod,EFermiR,EFermiL,
     &           kT,GFeta,NVOLT,contmp,WGFtmp,
     &           0,0,0)
           do i=1,NVOLT
              jc=jc+1
              contour(jc)=contmp(i)
              WGF(jc)=WGFtmp(i)
           end do              !i
           call memory('D','Z',2*nvolt,'setupcontour') 
           deallocate(contmp)
           deallocate(wgftmp)

        end if                 !finite voltage


        NEn=jc

        if(Node.eq.0) 
     .   write(*,*) NEn,' energy points'


        if(Node.eq.0) then
           fname = paste(slabel,'.CONTOUR')
           call io_assign( ju )
           open( ju, file=fname, status='unknown' )
           do i=1,NEn
              write(ju,'(4F8.3,i4)') contour(i),WGF(i),i 
           end do
           call io_close( ju )
        end if               

      call timer('SUcontour',2)

c ===================================================================
      return
      end subroutine setupcontour
c ===================================================================

!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------



      END MODULE m_ts_contour
