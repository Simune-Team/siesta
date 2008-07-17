      module m_projected_DOS
      private
      public :: projected_DOS

      CONTAINS
      
      subroutine projected_DOS()
      use precision
      use sparse_matrices
      USE siesta_options
      use atomlist, only: indxuo, no_s, no_u, no_l
      use fdf
      use sys, only: die
      use Kpoint_grid
      use parallel, only: IOnode
      use m_eo
      use m_spin,   only: nspin
      use m_gamma
      implicit none

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      real(dp) :: factor
      integer  :: nhist ! Number of histogram intervals in projected DOS
      real(dp) :: e1    ! Lower bound of energy range
      real(dp) :: e2    ! Upper bound of energy range
      real(dp) :: sigma ! Energy width used to convolute partial DOS

!------------------------------------------------------------------------- BEGIN
! Compute the projected density of states
      do_pdos = fdf_block('ProjectedDensityOfStates',bfdf)
      if ((isolve .ne. 0) .and. do_pdos) then
        write(6,*) 'siesta: ERROR: PDOS implemented only with diagon'
        do_pdos = .false.
      endif

      if (do_pdos) then
! Find the desired energy range
        if (.not. fdf_bline(bfdf,pline))
     $    call die('projected_DOS: ERROR in ' //
     $             'ProjectedDensityOfStates block')
        if ((fdf_bnvalues(pline) .lt. 3) .or.
     $      (fdf_bnnames(pline) .ne. 1))
     $    call die("Wrong format in PDOS block")
        factor = fdf_convfac( fdf_bnames(pline,1), 'Ry' )
        e1 = fdf_bvalues(pline,1) * factor
        e2 = fdf_bvalues(pline,2) * factor
        sigma = fdf_bvalues(pline,3) * factor
        nhist = fdf_bintegers(pline,1)
        if (IOnode) then
          write(6,'(a)') 'siesta: PDOS info: '
          write(6,'(a,3(f8.2,a),2x,i5)')
     $         'siesta: e1, e2, sigma, nhist: ',
     $         e1/eV,' eV',e2/eV,' eV',sigma/eV,' eV', nhist
        endif

        call pdos( no_s, nspin, nspin, no_l, maxnh,
     .             no_u, numh, listhptr, listh, H, S,
     .             e1, e2, sigma, nhist,
     .             gamma, xijo, indxuo, nkpnt, kpoint, kweight, eo,
     .             no_u)

      endif                     ! PDOS calc (do_pdos)
!--------------------------------------------------------------------------- END

      end subroutine projected_DOS
      end module m_projected_DOS
