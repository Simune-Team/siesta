      MODULE m_mixer
      private
      public :: mixer

      CONTAINS
 
      subroutine mixer( iscf )
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: iaorb, iphorb, lasto, no_u, no_l, indxuo ! CC RC
      use m_broyden_mixing
      use m_fire_mixing
      use m_spin,   only: nspin, h_spin_dim, SpOrb
      use m_pulay
      use parallel, only: IONode
#ifdef TRANSIESTA
      use m_ts_global_vars, only : TSiscf, TSrun
#endif

      implicit none

      real(dp), dimension(:,:), pointer :: Xin => null()
      real(dp), dimension(:,:), pointer :: Xout => null()  

      integer,  intent(in)  :: iscf

      real(dp)              :: dummy_dmax   ! Now local variable
      integer               :: iiscf
      logical               :: mmix    ! Same as mix. Used in pulayx

      external ::  mulliken

!----------------------------------------------------------- BEGIN

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
#ifdef TRANSIESTA
        if ( TSrun ) then
           iiscf=TSiscf
           if (maxsav .le. 0) then
             iiscf = 1
             if (TSiscf .ne. 1) mmix = .true.
           endif
        end if
#endif /* TRANSIESTA */

! CC RC  Changed in the following nspin by h_spin_dim
        if (fire_mix) then
          call fire_mixing( iscf, mix, no_l, maxnh,
     &                      numh(1:no_l), listhptr(1:no_l), h_spin_dim,
     $                      wmix, nkick, wmixkick, Xout, Xin,dummy_dmax)
        else if (broyden_maxit == 0) then
               call pulayx( iiscf, mmix, no_l, maxnh,
     &                 numh, listhptr, h_spin_dim, maxsav, wmix, nkick, 
     .                 wmixkick, XOut, XIn, dummy_dmax)
        else
          call broyden_mixing( iscf, mix, no_l, maxnh,
     &                         numh(1:no_l), listhptr(1:no_l),
     &                         h_spin_dim, wmix, nkick, wmixkick,
     $                         XOut, XIn, dummy_dmax)
        endif
      endif

!     Print populations at each SCF step, if requested
      ! Note that this is after mixing, which is not
      ! entirely correct. It should be moved to the top,
      ! or done somewhere else.

      if (muldeb) then 
         if (IONode)
     &        write (6,"(/a)")
     &        'siesta: Mulliken populations after mixing'
CC RC  osSO 
         if (SpOrb) then
          call moments( mullipop, na_u, no_u, maxnh, numh, listhptr,
     .                      listh, S, Dscf, isa, lasto, iaorb, iphorb,
     .                      indxuo )
         else
          call mulliken( mullipop, h_spin_dim, na_u, no_u, maxnh,
     &                 numh, listhptr, listh, S, Dscf, isa,
     &                 lasto, iaorb, iphorb )
         endif
CC RC 
      endif

      call timer( 'MIXER', 2 )

!-------------------------------------------------------- END
      END subroutine mixer

      End MODULE m_mixer


