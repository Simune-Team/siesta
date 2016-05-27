! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module tsread2
      interface
C FDN kscell and kdispl added as dummy 
      subroutine TSiohs2(task,finname, gamma, nua, nuotot, notot,
     .            nspin,maxnh,numh, listhptr, listh, H, S, qtot, temp,
     .               xij, indxuo, efs, cell, isa, lasto, xa,
     .            kscell, kdispl  )

      character         task*(*),finname*33
      logical           gamma
      integer           maxnh, nuotot, nspin, nua , notot
      double precision  qtot, temp, cell(3,3)
      integer, dimension (:), pointer:: listh, listhptr,
     .                        numh,indxuo,isa,lasto
      double precision, dimension (:,:), pointer:: H,xij,xa
      double precision, dimension (:), pointer:: S,efs
      integer kscell(3,3)
      real*8 kdispl(3)

      end subroutine TSiohs2

      subroutine TSiohs(task,finname, gamma, nua, nuotot, notot,
     .            nspin,maxnh,numh, listhptr, listh, H, S, qtot, temp,
     .               xij, indxuo, efs, cell, isa, lasto, xa,
     .            kscell, kdispl  )

      character         task*(*),finname*33
      logical           gamma
      integer           maxnh, nuotot, nspin, nua , notot
      double precision  qtot, temp, cell(3,3)
      integer, dimension (:), pointer:: listh, listhptr,
     .                        numh,indxuo,isa,lasto
      double precision, dimension (:,:), pointer:: H,xij,xa
      double precision, dimension (:), pointer:: S,efs
      integer kscell(3,3)
      real*8 kdispl(3)

      end subroutine TSiohs 

      end interface
      end module tsread2

