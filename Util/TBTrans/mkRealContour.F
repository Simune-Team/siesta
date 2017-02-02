! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
C ##################################################################
C ##      Generate (close to) real axis energy contour            ##
C ##           Transmission and DOS calculation                   ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################


      subroutine mkRealContour(joutfile,E1,E2,
     &     GFeta,nencont,zcontour,wgzcontour)

      use parallel, only : IOnode 

      implicit none     

c     PARAMETERS
      real*8 eV
      parameter ( eV = 1.d0 / 13.60580d0 )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT 
      integer joutfile          !OUT-FILE    
      integer nencont           ! No. contour points
      real*8 E1,E2              ! energy parameters 
      real*8 GFeta              ! state broadening in Ry


c OUTPUT
      complex*16 zcontour(nencont) ! points for GF
      complex*16 wgzcontour(nencont) ! weights on GF
      
      real*8 delta

      integer ic

c=======================================================================
c BEGIN:
c=======================================================================

 666  format(a10,f12.5,1X,f12.5,2X,f15.9,i4)

c
c     Simple Line
c
    
      delta = (E2 - E1) / nencont
      if (IOnode) then
        write(joutfile,*) 'contour:  Simple Line Contour:'
      end if
        
      do ic=1,nencont
         zcontour(ic) = dcmplx(E1+(ic-0.5d0)*delta,GFeta)
         wgzcontour(ic) = dcmplx(delta,0d0)
         if (IOnode) then
           write(joutfile,666) 'contour: ',
     &        dreal(zcontour(ic))/eV,GFeta/eV,delta/eV,ic
         end if ! IOnode
         
      end do                    !ic

C =========================================================
      return
      end
C =========================================================

