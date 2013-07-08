      MODULE m_compute_dm
      private
      public :: compute_dm

      CONTAINS

      subroutine compute_dm( iscf )
      use precision
      use units, only: eV
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: qa, lasto, iphorb, iaorb, no_u, no_s, indxuo, 
     &                    qtot, no_l
      use sys, only: die
      use Kpoint_grid
      use m_energies,   only: Ebs, DEHarr, Ecorrec, Ef, Entropy
      use m_rmaxh
      use m_eo
      use m_gamma
      use parallel, only: IONode
#ifdef MPI
      use m_mpi_utils,     only:  globalize_sum
#endif
#ifdef CDF
      use iodmhs_netcdf, only: write_dmh_netcdf
#endif
      use m_ordern,        only : ordern
      use m_steps,         only : istp
      use m_spin,          only : nspin, qs, efs
#ifdef TRANSIESTA
      use m_ts_global_vars, only : TSiscf, TSinit, TSrun
      use m_ts_method
      use m_ts_options,     only : TSmode, ElLeft, ElRight
      use m_ts_options,     only : no_BufL, no_BufR
      use m_ts_full,        only : transiesta_full
      use m_ts_tri,         only : ts_tri_init
      use m_ts_tri,         only : transiesta_tri
      use m_ts_sparse,      only : ts_sparse_init
      use m_ts_charge,      only : ts_print_charges
      use m_ts_charge,      only : TS_RHOCORR_METHOD
      use m_ts_charge,      only : ts_charge_correct
#endif /* TRANSIESTA */
      implicit none
!     Input variables
      integer, intent(in) :: iscf
!     Local variables
      integer             :: io, ispin, nnz
      real(dp)            :: const ! Temporary factor used in
                             ! calculation of Harris forces. e1>e2 to
                             ! signal that we do not want DOS weights
      real(dp), parameter :: e1 = 1.0_dp, e2 = -1.0_dp
#ifdef MPI
      real(dp)            :: buffer1
#endif

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE compute_dm' )
#endif
      call timer( 'compute_dm', 1 )

! Save present DM here, instead of in setup_hamiltonian...

      do ispin = 1,nspin
         do io = 1,maxnh
            Dold(io,ispin) = Dscf(io,ispin)
         enddo
      enddo

!     Solve eigenvalue problem 
      if (isolve .eq. SOLVE_DIAGON) then
        call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh, 
     &              H, S, qtot, fixspin, qs, temp, e1, e2,
     &              gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &              occtol, iscf, neigwanted)
        Ecorrec = 0.0_dp
      elseif (isolve .eq. SOLVE_ORDERN) then
        if (.not. gamma) call die("Cannot do O(N) with k-points.")
        call ordern(usesavelwf, ioptlwf, na_u, no_u, no_l, lasto,
     &               isa, qa, rcoor, rmaxh, ucell, xa, iscf,
     &               istp, ncgmax, etol, eta, qtot, maxnh, numh,
     &               listhptr, listh, H, S, chebef, noeta, rcoorcp,
     &               beta, pmax, Dscf, Escf, Ecorrec, nspin, qs )
        Entropy = 0.0_dp
#ifdef TRANSIESTA
! TSS Begin
      elseif (TSmode .and. TSinit) then
        call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh,
     &              H, S, qtot, fixspin, qs, temp, e1, e2,
     &              gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &              occtol, iscf, neigwanted)

        Ecorrec = 0._dp

      else if (TSrun) then
        call timer('TS',1)

#ifdef TRANSIESTA_DEBUG
        ! When we do a transiesta DEBUG run, we will write out
        ! all matrices: SigmaL, SigmaR, Gf^-1, Gf
        ! To compare we "force" both the tri-diagonal and the
        ! regular method
        if ( TSiscf == 1 ) then
           ! This is a first Transiesta run, we need to initialize
           ! the sparsity patterns...
           call ts_sparse_init(Gamma, block_dist, sparse_pattern,
     &             na_u, lasto)
           call ts_tri_init()
           call ts_print_charges(ElLeft, ElRight, no_BufL, no_BufR,
     &          block_dist, sparse_pattern,nspin, maxnh, Dscf, S)
        end if
        call transiesta_full(nspin,Gamma, block_dist,
     &       sparse_pattern,ucell, no_u, na_u, lasto, xa, maxnh,
     &       xijo, H, S, Dscf,Escf, Ef, TSiscf, Qtot)
        call transiesta_tri(nspin,Gamma, block_dist,
     &       sparse_pattern,ucell, no_u, na_u, lasto, xa, maxnh,
     &       xijo, H, S, Dscf,Escf, Ef, TSiscf, Qtot)
#else
        if ( TSiscf == 1 ) then
           ! This is a first Transiesta run, we need to initialize
           ! the sparsity patterns...
           call ts_sparse_init(Gamma, block_dist, sparse_pattern,
     &             na_u, lasto)
           if ( ts_method == TS_SPARSITY_TRI ) then
              call ts_tri_init()
           end if
           call ts_print_charges(ElLeft, ElRight, no_BufL, no_BufR,
     &          block_dist, sparse_pattern,nspin, maxnh, Dscf, S)
        end if
        if ( ts_method == TS_SPARSITY ) then
           call transiesta_full(nspin,Gamma, block_dist,
     &             sparse_pattern,ucell, no_u, na_u, lasto, xa, maxnh,
     &          xijo, H, S, Dscf,Escf, Ef, TSiscf, Qtot)
        else if ( ts_method == TS_SPARSITY_TRI ) then
           call transiesta_tri(nspin,Gamma, block_dist,
     &          sparse_pattern,ucell, no_u, na_u, lasto, xa, maxnh,
     &          xijo, H, S, Dscf,Escf, Ef, TSiscf, Qtot)
        end if
#endif

        ! We do the charge correction of the transiesta
        ! computation here (notice that the routine will automatically
        ! return if no charge-correction is requested)
        call ts_charge_correct(ElLeft, ElRight,
     &       no_BufL, no_BufR, block_dist, sparse_pattern, 
     &       nspin, maxnh, Dscf, Escf, S, Qtot, TS_RHOCORR_METHOD)

           
        call timer('TS',2)

        if (TSiscf .eq. 1) call timer( 'TS', 3 ) 
        if (TSiscf .eq. 1) call timer( 'TS_calc', 3 ) 

        Ecorrec = 0._dp
        Entropy = 0.0_dp
! TSS End
#endif /* TRANSIESTA */
      else
        call die('siesta: ERROR: wrong solution method')
      endif

! ..................
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
     &    'siesta: Mulliken populations before mixing'
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

      call timer( 'compute_dm', 2 )
#ifdef DEBUG
      call write_debug( '  POS compute_dm' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine compute_dm
      END MODULE m_compute_dm
