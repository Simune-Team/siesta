      MODULE m_mixer
      private
      public :: mixer

      CONTAINS
 
      subroutine mixer( iscf, dDmax )
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: iaorb, iphorb, lasto, no_u, no_l
      use m_broyden_mixing
      use m_fire_mixing
      use m_steps
      use m_spin,   only: nspin
      use m_pulay
      use m_iodm,   only: write_dm
      use parallel, only: IONode
#ifdef MPI
      use m_mpi_utils, only: globalize_max
#endif
#ifdef CDF
      use iodm_netcdf, only: write_dm_netcdf
#endif

      implicit none

      integer,  intent(in)   :: iscf
      real(dp), intent(out) :: dDmax     ! Max. change in density matrix elements

      integer :: iiscf
      logical :: mmix       ! Same as mix. Used in pulayx

#ifdef MPI
      real(dp):: buffer1   ! Temporary for the result of a globalize operation
#endif
      external ::  iodm, mulliken

!------------------------------------------------------------------------- BEGIN
! Mix input and output energy-density and density matrices ............
! Following line for using and saving the density matrix without mix ..
          if (wmix.ne.0._dp) then
! Pulay or Broyden mixing
            mmix  = mix
            iiscf = iscf
            if (maxsav .le. 0) then
              iiscf = 1
              if (iscf .ne. 1) mmix = .true.
            endif
            if (fire_mix) then
               call fire_mixing(iscf, mix, no_l, maxnh,
     .                 numh(1:no_l), listhptr(1:no_l), nspin,
     $                 wmix, nkick, wmixkick, Dscf, Dold, dDmax)
            else if (broyden_maxit == 0) then
               call pulayx( pulfile, iiscf, mmix, no_l, maxnh,
     .                   numh, listhptr, nspin, maxsav, wmix, nkick, 
     .                   wmixkick, Dscf, Dold, dDmax)
            else
               call broyden_mixing(iscf, mix, no_l, maxnh,
     .                 numh(1:no_l), listhptr(1:no_l), nspin,
     $                 wmix, nkick, wmixkick, Dscf, Dold, dDmax)
            endif
          endif

! Ensure that dDmax is the same on all nodes for convergence test/output
#ifdef MPI
          call globalize_max(dDmax,buffer1)
          dDmax = buffer1
#endif
! ...................

! Print populations at each SCF step, if requested, after mixing ......

          if (muldeb) then 
             write (6,"(/a)")
     .             'siesta: Mulliken populations after mixing'
             call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .                      numh, listhptr, listh, S, Dscf, isa,
     .                      lasto, iaorb, iphorb )
          endif
! ..................

! Save density matrix on disk, after mixing ...........................
#ifdef CDF
          if (writedm_cdf_history) then
             call write_dm_netcdf( no_l, maxnh, nspin, Dscf,
     $                             overwrite=.false.)
          else if (writedm_cdf) then
             call write_dm_netcdf( no_l, maxnh, nspin, Dscf,
     $                             overwrite=.true. )
          endif
#endif
    
          if (writedm) then
            if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
              if (istp.eq.1)
!!     .        call iodm( 'write', maxnh, no_l, nspin,
!!     .                   numh, listhptr, listh, Dscf, found )
     .        call write_dm (maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf)
            else
               call write_dm (maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf)
            endif
          endif
!--------------------------------------------------------------------------- END
      END subroutine mixer

      End MODULE m_mixer


