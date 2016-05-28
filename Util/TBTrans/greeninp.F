! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module greeninp
      interface 
      subroutine read_green(joutfile,jgfu,
     .     EFermi,NEn,contour,wgf,
     .     nua,NA1,NA2,nq,ng,wq,qi,kxy,iq)

      integer joutfile          ! out file
      integer jgfu               !unit of gf-file
      real*8 efermi             !The required Fermi energy
      integer NEn,ng,iq
      complex*16 contour(NEn),wgf(NEn)
      integer nua,NA1,NA2       !no. atoms in uc and no. repetitions in A1,A2
      integer nq                !no. q-points
      real*8, dimension (:,:), pointer:: qi ! the shifted q-points
      real*8 kxy(2) ! the actual extended cell k point.
      real*8, dimension (:), pointer:: wq

      end subroutine read_green

      end interface
      end module greeninp

