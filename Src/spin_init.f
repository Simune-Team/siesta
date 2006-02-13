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
      subroutine spin_init(nspin)

C integer nspin            : Spin polarization
C real*8 cell(3,3)         : (Super) lattice vectors CELL(ixyz,ivector)
C                            (in Bohr)
C integer ncells           : Number of unit cells in supercell
C **********************************************************************

C
C  Modules
C
      use precision
      use fdf

      implicit none

      integer, intent(out)  :: nspin

      logical  noncol, sppol

      nspin = 1

      sppol  = fdf_boolean('SpinPolarized',.false.)
      noncol = fdf_boolean('NonCollinearSpin',.false.)

      if (noncol) then
        nspin = 4
      elseif (sppol) then
        nspin = 2
      else 
        nspin = 1
      endif

      write(6,'(a,4x,i1)') 
     . 'redata: Number of spin components        = ',nspin

      end
