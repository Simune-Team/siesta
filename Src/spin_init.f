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

      sppol  = fdf_get('SpinPolarized',.false.)
      noncol = fdf_get('NonCollinearSpin',.false.)

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
