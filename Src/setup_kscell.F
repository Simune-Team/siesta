! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine setup_kscell( cell, kscell, displ, firm_displ )

c ***************** INPUT **********************************************
c real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
c ***************** OUTPUT *********************************************
c integer kscell(3,3)  : Supercell reciprocal of k-grid unit cell
c real*8  displ(3)     : Grid origin in k-grid-vector coordinates:
c logical firm_displ   : User-specified displacements

c   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
c   If both are present, kgrid_Monkhorst_Pack has priority. If none is
c   present, the cutoff default is zero, producing only the gamma point.
c   Examples of fdf data specifications:
c     kgrid_cutoff  50. Bohr
c     %block kgrid_Monkhorst_Pack  # Defines kscell and displ
c     4  0  0   0.50               # (kscell(i,1),i=1,3), displ(1)
c     0  4  0   0.50               # (kscell(i,2),i=1,3), displ(2)
c     0  0  4   0.50               # (kscell(i,3),i=1,3), displ(3)
c     %endblock kgrid_Monkhorst_Pack
c **********************************************************************
C
C  Modules
C
      use precision,  only : dp
      use units,      only : Ang
      use parallel,   only : Node
      use m_minvec,   only : minvec
      use fdf

      implicit          none

C Passed variables
      real(dp), intent(in)   :: cell(3,3)
      integer, intent(out)   :: kscell(3,3)
      real(dp), intent(out)  :: displ(3)
      logical, intent(out)   :: firm_displ

C Internal variables
      integer           i, j,  factor(3,3), expansion_factor
      real(dp)          scmin(3,3),  vmod
      real(dp)          ctransf(3,3)
      logical           genlogic

      real(dp), parameter :: defcut = 0.0_dp
      integer, dimension(3,3), parameter :: unit_matrix = 
     .                    reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/))

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

#ifdef DEBUG
      call write_debug( '    PRE setup_kscell' )
#endif

      genlogic = fdf_block('kgrid_Monkhorst_Pack',bfdf)
      if ( genlogic ) then
        do i = 1,3
          fdf_bline(bfdf,pline)
          kscell(1,i) = fdf_bintegers(pline,1)
          kscell(2,i) = fdf_bintegers(pline,2)
          kscell(3,i) = fdf_bintegers(pline,3)
          if ( fdf_bnvalues(pline) > 3 ) then
             displ(i) = fdf_bvalues(pline,4)
          else
             displ(i) = 0._dp
          end if
        enddo
        firm_displ = .true.

      else

        cutoff = fdf_physical('kgrid_cutoff',defcut,'Bohr')
        if (Node .eq. 0) then
          write(6,"(a,f10.5)") "Kgrid cutoff input: ", cutoff
        endif
        firm_displ = .false.
        
        ! Find equivalent rounded unit-cell
        call minvec( cell, scmin, ctransf )

        expansion_factor = 1
        do j = 1,3
           factor(j,1:3) = 0
           vmod = sqrt(dot_product(scmin(1:3,j),scmin(1:3,j)))
           factor(j,j) = int(2.0_dp*cutoff/vmod) + 1
           if (mod(factor(j,j),2) == 0) then
              displ(j) = 0.5_dp
           else
              displ(j) = 0.0_dp
           endif
           expansion_factor = expansion_factor * factor(j,j)
        enddo
        ! Generate actual supercell skeleton
        kscell = matmul(ctransf, factor)
        ! Avoid confusing permutations
        if (expansion_factor == 1) then
           kscell = unit_matrix
        endif
      endif

#ifdef DEBUG
      call write_debug( '    POS setup_kscell' )
#endif

      end subroutine setup_kscell

