      MODULE m_compute_dm
      private
      public :: compute_dm

      logical, public, save :: PreviousCallDiagon=.false.

      CONTAINS

      subroutine compute_dm( iscf )
      use precision
      use units, only: eV
      USE siesta_options
      use class_dSpData1D, only : val
      use sparse_matrices
      use siesta_geom
      use atomlist, only: qa, lasto, iphorb, iaorb, no_u, no_s, indxuo, 
     &                    qtot, no_l
      use sys, only: die
      use Kpoint_grid
      use m_energies,   only: Ebs, Ecorrec, Ef, Entropy
      use m_rmaxh
      use m_eo
      use m_gamma
      use parallel, only: IONode
#ifdef CDF
      use iodmhs_netcdf, only: write_dmh_netcdf
#endif
      use m_dminim,        only : dminim
      use m_zminim,        only : zminim
      use m_ordern,        only : ordern
      use m_steps,         only : istp
      use m_spin,          only : nspin, qs, efs
      use m_normalize_dm,  only : normalize_dm
#ifdef TRANSIESTA
      use m_ts_global_vars, only : TSmode, TSinit, TSrun
      use m_transiesta,     only : transiesta
#endif /* TRANSIESTA */
      implicit none
!     Input variables
      integer, intent(in) :: iscf
!     Local variables
      logical             :: CallDiagon
      integer             :: nnz
      real(dp), pointer :: H_kin(:)

      ! e1>e2 to signal that we do not want DOS weights
      real(dp), parameter :: e1 = 1.0_dp, e2 = -1.0_dp

!-------------------------------------------------------------- BEGIN
      call timer( 'compute_dm', 1 )

! Save present DM here, instead of in setup_hamiltonian...

!$OMP parallel workshare default(shared)
      Dold(:,:) = Dscf(:,:)
!$OMP end parallel workshare

!     Solve eigenvalue problem 

      ! Here we decide if we want to calculate one or more SCF steps by
      ! diagonalization before proceeding with the OMM routine
      CallDiagon=.false.
      if (isolve .eq. SOLVE_MINIM) then
        if (istp .eq. 1) then
          if ((iscf .le. call_diagon_first_step) .or.
     &        (call_diagon_first_step<0)) CallDiagon=.true.
        else
          if ((iscf .le. call_diagon_default) .or.
     &        (call_diagon_default<0)) CallDiagon=.true.
        endif
      endif

      if ((isolve .eq. SOLVE_DIAGON) .or. (CallDiagon)) then
        call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh, 
     &              H, S, qtot, fixspin, qs, temp, e1, e2,
     &              gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &              occtol, iscf, neigwanted)
        Ecorrec = 0.0_dp
        PreviousCallDiagon=.true.
      elseif (isolve .eq. SOLVE_ORDERN) then
        if (.not. gamma) call die("Cannot do O(N) with k-points.")
        call ordern(usesavelwf, ioptlwf, na_u, no_u, no_l, lasto,
     &               isa, qa, rcoor, rmaxh, ucell, xa, iscf,
     &               istp, ncgmax, etol, eta, qtot, maxnh, numh,
     &               listhptr, listh, H, S, chebef, noeta, rcoorcp,
     &               beta, pmax, Dscf, Escf, Ecorrec, nspin, qs )
        Entropy = 0.0_dp
      elseif (isolve .eq. SOLVE_MINIM) then
        if (nspin==4) call die('ERROR: Non-collinear spin calculations
     &                          not yet implemented with the OMM!')
        H_kin => val(H_kin_1D)
        if (gamma) then
          call dminim(.false., PreviousCallDiagon, iscf, istp, no_l,
     &                 nspin, no_u, maxnh, numh, listhptr, listh, Dscf,
     &                 eta, qs, H, S, H_kin)
        else
          call zminim(.false., PreviousCallDiagon, iscf, istp, no_l,
     &                 nspin, no_u, maxnh, numh, listhptr, listh, Dscf,
     &                 eta, qs, no_s, xijo, indxuo, nkpnt, kpoint,
     &                 kweight, H, S, H_kin)
        end if
        Ecorrec = 0.0_dp
        Entropy = 0.0_dp
        PreviousCallDiagon=.false.
#ifdef TRANSIESTA
      elseif (TSmode .and. TSinit) then
        call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh,
     &              H, S, qtot, fixspin, qs, temp, e1, e2,
     &              gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &              occtol, iscf, neigwanted)

        Ecorrec = 0._dp

      else if (TSrun) then

         call transiesta(iscf,nspin, block_dist, sparse_pattern,
     &        Gamma, ucell, nsc, isc_off, no_u, na_u, lasto, xa, maxnh,
     &        H, S, Dscf, Escf, Ef, Qtot, .false. )

         Ecorrec = 0._dp
         Entropy = 0.0_dp
#endif /* TRANSIESTA */
      else
        call die('siesta: ERROR: wrong solution method')
      endif

#ifdef CDF
      if ( writedmhs_cdf_history) then
        call write_dmh_netcdf( no_l, maxnh, nspin, Dold, H, Dscf )
      else if (writedmhs_cdf) then
        call write_dmh_netcdf( no_l, maxnh, nspin, Dold, H, Dscf,
     &                         overwrite=.true. )
      endif
#endif

!     Print populations at each SCF step if requested before mixing ......
      if (muldeb) then
        if (ionode) write (6,"(/a)")
     &    'siesta: Mulliken populations (DM_out)'
        call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     &                 numh, listhptr, listh, S, Dscf, isa,
     &                 lasto, iaorb, iphorb )
      endif

! Write orbital indexes. JMS Dec.2009

#ifdef TRANSIESTA
      if (IOnode .and. iscf==1.and..not. TSrun) then
         call write_orb_indx( na_u, na_s, no_u, no_s, isa, xa,
     .        iaorb, iphorb, indxuo, nsc, ucell )
      endif
#else
      if (IOnode .and. iscf==1) then
         call write_orb_indx( na_u, na_s, no_u, no_s, isa, xa,
     .        iaorb, iphorb, indxuo, nsc, ucell )
      endif
#endif

!     Normalize density matrix to exact charge
!     Placed here for now to avoid disturbing EHarris
#ifdef TRANSIESTA      
      if ( .not. TSrun ) call normalize_dm( first= .false. )
#else
      call normalize_dm(first=.false.)
#endif
      
      call timer( 'compute_dm', 2 )

!----------------------------------------------------------------- END
      END subroutine compute_dm
      END MODULE m_compute_dm
