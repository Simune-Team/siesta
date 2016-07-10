! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      

      module hamil2
      interface 
      subroutine setupkham2(joutfile,tinit,kpoint,ispin,hsfile,
     &  lasto,Hk,Sk,nua,nuo,nspin) 
      integer joutfile        
      logical tinit
      real*8 kpoint(3) 
      integer ispin            
      character*33 hsfile    
      integer, dimension (:), pointer:: lasto
      complex*16, dimension (:), pointer:: Hk, Sk
      integer nua,nuo,nspin
      end subroutine setupkham2

      end interface
      end module hamil2

