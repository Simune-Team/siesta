! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_projected_DOS

      use precision
      
      implicit none

      private

      public :: init_projected_DOS, projected_DOS

      CONTAINS
      
      subroutine init_projected_DOS()

      USE siesta_options
      use fdf,         only: fdf_block, block_fdf
      ! This is to get the reference kpoints in case PDOS.kgrid* has not
      ! been specified
      use kpoint_pdos_m, only: setup_kpoint_pdos
      use parallel,    only: IOnode, Nodes
      use siesta_geom, only: ucell
      use m_spin,      only: spin

      type(block_fdf) :: bfdf
      
!-------------------------------------------------------------------------BEGIN
!     Compute the projected density of states
      do_pdos = fdf_block('ProjectedDensityOfStates', bfdf)
      if ( .not. do_pdos ) return
         
      if (isolve.ne.SOLVE_DIAGON) then
         if (.not.((isolve.eq.SOLVE_MINIM).and.
     .        minim_calc_eigenvalues)) then
            if (IONode) then
               write(6,*)
     .              'siesta: ERROR: PDOS implemented only with diagon'
            endif
            do_pdos = .false.
         endif
      endif

      if ( .not. do_pdos ) return

!     Determine whether the projected density of states is to be computed
!     on a different grid to the SCF calculation
      call setup_kpoint_pdos( ucell )
      
!---------------------------------------------------------------------------END
      
      end subroutine init_projected_DOS

      subroutine projected_DOS()

      use sparse_matrices
      USE siesta_options
      use alloc,       only : re_alloc
      use atomlist,    only : indxuo, no_s, no_u, no_l
      use fdf
      use sys,         only : die
      use kpoint_scf_m, only: kpoints_scf
      use Kpoint_pdos_m, only: kpoints_pdos, gamma_pdos
      use parallel,    only: IOnode
      use m_eo
      use m_spin,      only: h_spin_dim, spinor_dim
      use m_gamma
      use units,       only: eV
      
      implicit none

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      real(dp) :: factor
      logical  :: dummy ! Logical to hold return value from call to fdf_block
      integer  :: nhist ! Number of histogram intervals in projected DOS
      real(dp) :: e1    ! Lower bound of energy range
      real(dp) :: e2    ! Upper bound of energy range
      real(dp) :: sigma ! Energy width used to convolute partial DOS

!------------------------------------------------------------------------- BEGIN
! Compute the projected density of states

      if ( .not. do_PDOS ) return

      ! Call fdf_block to get iu - presence has already been tested in init_projected_PDOS
      dummy = fdf_block('ProjectedDensityOfStates',bfdf)
! Find the desired energy range
      if (.not. fdf_bline(bfdf,pline))
     $     call die('projected_DOS: ERROR in ' //
     $     'ProjectedDensityOfStates block')
      if ((fdf_bnreals(pline) .lt. 3) .or.
     $     (fdf_bnnames(pline) .ne. 1))
     $   call die("Wrong format in PDOS block, not enough reals/names")
      factor = fdf_convfac( fdf_bnames(pline,1), 'Ry' )
      e1 = fdf_breals(pline,1) * factor
      e2 = fdf_breals(pline,2) * factor
      sigma = fdf_breals(pline,3) * factor
      nhist = fdf_bintegers(pline,1)
      if (IOnode) then
         write(6,'(a)') 'siesta: PDOS info: '
         write(6,'(a,3(f8.2,a),2x,i5)')
     $        'siesta: e1, e2, sigma, nhist: ',
     $        e1/eV,' eV',e2/eV,' eV',sigma/eV,' eV', nhist
       end if

      ! If the k points have been set specifically for the PDOS then use this set
      if ( kpoints_pdos%N > kpoints_scf%N ) then
        call re_alloc(eo,1,no_u,1,spinor_dim,1,kpoints_pdos%N,name="eo",
     .      routine="projected_dos")
      end if
      
      call pdos( no_s, h_spin_dim, spinor_dim, no_l, 
     .    maxnh, no_u, numh, listhptr, listh, H, S,
     .    e1, e2, sigma, nhist, gamma_pdos, xijo, indxuo, 
     .    kpoints_pdos%N, kpoints_pdos%k, kpoints_pdos%w, eo,
     .    no_u)

      end subroutine projected_DOS

      end module m_projected_DOS
