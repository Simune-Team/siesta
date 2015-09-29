      MODULE m_mixer
      private
      public :: mixer

      CONTAINS
 
      subroutine mixer( iscf )

      use precision,        only: dp
      use siesta_options,   only: mixH
      use siesta_options,   only: mullipop, muldeb
      use sparse_matrices,  only: Dold, Dscf, Hold, H, S
      use sparse_matrices,  only: maxnh, numh, listhptr, listh
      use siesta_geom,      only: na_u, isa
      use atomlist,         only: iaorb, iphorb, lasto, no_u, no_l

      use m_mixing_scf, only: scf_mix
      use m_mixing, only: mixing
      
#ifdef OLD_MIXING
      use m_broyden_mixing, only: broyden_mixing
      use m_fire_mixing,    only: fire_mixing
#endif

      use m_spin,           only: nspin
      use parallel,         only: IONode

      implicit none

      real(dp), pointer :: Xin(:,:), Xout(:,:)

      integer,  intent(in)  :: iscf

      real(dp)              :: dummy_dmax   ! Now local variable
      integer               :: iiscf
      logical               :: mmix    ! Same as mix. Used in pulayx

      external :: mulliken

!-------------------------------------------------------------------- BEGIN

      call timer( 'MIXER', 1 )

      if ( mixH ) then
         ! Mix Hamiltonian
         Xin => Hold
         Xout => H
      else
         ! Mix density matrix
         Xin => Dold
         Xout => Dscf
      end if

#ifdef OLD_MIXING
      ! Convoluted legacy logic for the Pulay mixer
      ! iiscf is the "effective SCF step"
      ! To be cleaned up
      
      mmix  = mix_first_scf_step
      iiscf = iscf
      if (maxsav .le. 0) then
         iiscf = 1
         if (iscf .ne. 1) mmix = .true.
      end if
      
      if (fire_mix) then
         call fire_mixing( iscf, mix_first_scf_step, no_l, maxnh,
     &        numh(1:no_l), listhptr(1:no_l), nspin,
     $        wmix, nkick, wmixkick, Xout, Xin,dummy_dmax)
      else if (broyden_maxit == 0) then
         call pulayx( iiscf, mmix, no_l, maxnh,
     &        numh, listhptr, nspin, maxsav, wmix, nkick, 
     &        wmixkick, XOut, XIn, dummy_dmax)
      else
         call broyden_mixing( iscf, mix_first_scf_step, no_l, maxnh,
     &        numh(1:no_l), listhptr(1:no_l),
     &        nspin, wmix, nkick, wmixkick,
     $        XOut, XIn, dummy_dmax)
      endif

#endif 

      ! Call mixing routine
      call mixing( scf_mix, iscf, maxnh, nspin, Xin, Xout )

      ! Print populations at each SCF step, if requested
      ! Note that this is after mixing, which is not
      ! entirely correct. It should be moved to the top,
      ! or done somewhere else.

      if (muldeb) then 
         if (IONode)
     &        write (6,"(/a)")
     &        'siesta: Mulliken populations after mixing'
        call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     &                 numh, listhptr, listh, S, Dscf, isa,
     &                 lasto, iaorb, iphorb )
      endif

      call timer( 'MIXER', 2 )

!-------------------------------------------------------- END
      END subroutine mixer

      End MODULE m_mixer


