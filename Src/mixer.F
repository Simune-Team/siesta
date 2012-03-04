      MODULE m_mixer
      private
      public :: mixer

      CONTAINS
 
      subroutine mixer( iscf )

      use precision,        only: dp
      use siesta_options,   only: mix_first_scf_step, mixH, wmix
      use siesta_options,   only: wmixkick, fire_mix, nkick
      use siesta_options,   only: maxsav, broyden_maxit
      use siesta_options,   only: mullipop, muldeb
      use sparse_matrices,  only: Dold, Dscf, Hold, H, S
      use sparse_matrices,  only: maxnh, numh, listhptr, listh
      use siesta_geom,      only: na_u, isa
      use atomlist,         only: iaorb, iphorb, lasto, no_u, no_l
      use m_broyden_mixing, only: broyden_mixing
      use m_fire_mixing,    only: fire_mixing
      use m_pulay,          only: pulayx
      use m_spin,           only: nspin
      use parallel,         only: IONode
      use m_save_density_matrix,  only: save_density_matrix
#ifdef TRANSIESTA
      use m_energies,   only : ef
      use m_ts_global_vars, only: TSrun, TSiscf
      use sparse_matrices,  only: Escf
#endif /* TRANSIESTA */
      implicit none

      real(dp), dimension(:,:), pointer :: Xin => null()
      real(dp), dimension(:,:), pointer :: Xout => null()  

      integer,  intent(in)  :: iscf

      integer               :: iiscf
      logical               :: mmix    ! Same as mix. Used in pulayx
      real(dp)              :: dDmax   ! Dummy for calls

      external ::  mulliken

!-------------------------------------------------------------------- BEGIN

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

      call save_density_matrix("mixed")

      call timer( 'MIXER', 2 )
!------------------------------------------------------------------------ END
      END subroutine mixer

      End MODULE m_mixer


