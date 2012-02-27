      MODULE m_mixer
      private
      public :: mixer

      CONTAINS
 
      subroutine mixer( iscf )

      use precision,        only: dp
      use siesta_options,   only: mix_first_scf_step, mixH, wmix
      use siesta_options,   only: wmixkick, fire_mix, nkick
      use siesta_options,   only: maxsav, idyn, broyden_maxit
      use siesta_options,   only: mullipop, muldeb
      use siesta_options,   only: writedm, writedm_cdf
      use siesta_options,   only: writedm_cdf_history
      use sparse_matrices,  only: Dold, Dscf, Hold, H, S
      use sparse_matrices,  only: maxnh, numh, listhptr, listh
      use siesta_geom,      only: na_u, isa
      use atomlist,         only: iaorb, iphorb, lasto, no_u, no_l
      use m_broyden_mixing, only: broyden_mixing
      use m_fire_mixing,    only: fire_mixing
      use m_steps,          only: istp
      use m_spin,           only: nspin
      use m_pulay,          only: pulayx
      use m_iodm,           only: write_dm
      use parallel,         only: IONode
#ifdef CDF
      use iodm_netcdf, only: write_dm_netcdf
#endif
#ifdef TRANSIESTA
      use m_energies,   only : ef
      use m_ts_iodm, only : write_ts_dm
      use m_ts_global_vars, only: TSrun, TSiscf
#endif /* TRANSIESTA */
      implicit none

      real(dp), dimension(:,:), pointer :: Xin => null()
      real(dp), dimension(:,:), pointer :: Xout => null()  

      integer,  intent(in)  :: iscf

      integer               :: iiscf
      logical               :: mmix    ! Same as mix. Used in pulayx
      real(dp)              :: dDmax   ! Dummy for calls

      external ::  mulliken

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE mixer' )
#endif

      call timer( 'MIXER', 1 )

      if (mixH) then
         ! Mix Hamiltonian
         Xin => Hold
         Xout => H
      else
         ! Mix density matrix
         Xin => Dold
         Xout => Dscf
      endif


      if (wmix.ne.0._dp) then

        ! Convoluted legacy logic for the Pulay mixer
        ! iiscf is the "effective SCF step"
        ! To be cleaned up

        mmix  = mix_first_scf_step
        iiscf = iscf
        if (maxsav .le. 0) then
          iiscf = 1
          if (iscf .ne. 1) mmix = .true.
        endif
#ifdef TRANSIESTA
        if ( TSrun ) then
           iiscf=TSiscf
           if (maxsav .le. 0) then
             iiscf = 1
             if (TSiscf .ne. 1) mmix = .true.
           endif
        end if
#endif /* TRANSIESTA */

        if (fire_mix) then
          call fire_mixing( iscf, mix_first_scf_step, no_l, maxnh,
     &                      numh(1:no_l), listhptr(1:no_l), nspin,
     $                      wmix, nkick, wmixkick, Xout, Xin, dDmax)
        else if (broyden_maxit == 0) then
               call pulayx( iiscf, mmix, no_l, maxnh,
     &                 numh, listhptr, nspin, maxsav, wmix, nkick, 
     .                 wmixkick, XOut, XIn, dDmax)
        else
          call broyden_mixing( iscf, mix_first_scf_step, no_l, maxnh,
     &                         numh(1:no_l), listhptr(1:no_l),
     &                         nspin, wmix, nkick, wmixkick,
     $                         XOut, XIn, dDmax)
        endif
      endif

!     Print populations at each SCF step, if requested, after mixing ......
      if (muldeb) then 
        write (6,"(/a)")
     &    'siesta: Mulliken populations after mixing'
        call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     &                 numh, listhptr, listh, S, Dscf, isa,
     &                 lasto, iaorb, iphorb )
      endif

#ifdef CDF
!     Save density matrix on disk, after mixing, to
!     be used for re-starting the SCF cycle.
      if (writedm_cdf_history) then
         call write_dm_netcdf( no_l, maxnh, nspin, Dscf,
     &                         overwrite=.false. )
      else if (writedm_cdf) then
         call write_dm_netcdf( no_l, maxnh, nspin, Dscf,
     &                         overwrite=.true. )
      endif
#endif
    
#ifndef TRANSIESTA
      if (writedm) then
        if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
          if (istp.eq.1)
     &    call write_dm (maxnh, no_l, nspin,
     &               numh, listhptr, listh, Dscf)
        else
          call write_dm (maxnh, no_l, nspin,
     &               numh, listhptr, listh, Dscf)
        endif
      endif ! writedm
#else /* TRANSIESTA */
! TSS Begin
      if (writedm) then
        if (.not.TSrun) then   !TSS save DM
          if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
            if (istp.eq.1)
     &        call write_dm( maxnh, no_l, nspin, numh,
     &                       listhptr, listh, Dscf )
          else
            call write_dm( maxnh, no_l, nspin, numh,
     &                     listhptr, listh, Dscf)
          endif
        else
!         TSS write DME
             if ((idyn .eq. 6) .or. (idyn .eq. 7)) then
               if (istp .eq.1)
     .           call write_ts_dm (maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, Escf,ef)
             else
                call write_ts_dm (maxnh, no_l, nspin,
     .                   numh, listhptr, listh, Dscf, Escf,ef)
             endif            ! idyn.eq.6
        endif ! TSrun
      endif ! writedm
#endif /* TRANSIESTA */

      call timer( 'MIXER', 2 )
#ifdef DEBUG
      call write_debug( '  POS mixer' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine mixer

      End MODULE m_mixer


