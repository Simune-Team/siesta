! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module cont
      interface 
      subroutine setupcontour(joutfile,NEn,
     .     EFermi0,EFermiL,EFermiR,
     .     WGF,contour,ncontour,ncontour0,nvolt)

      integer joutfile,nen,ncontour,ncontour0,nvolt
      real*8 EFermi0,EFermiL,EFermiR  
      complex*16, dimension (:), pointer:: contour,wgf
      end subroutine setupcontour

      end interface
      end module cont

