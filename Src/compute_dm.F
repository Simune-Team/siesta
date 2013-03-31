      MODULE m_compute_dm
      private
      public :: compute_dm

      logical, public, save :: PreviousCallDiagon=.false.

      CONTAINS

      subroutine compute_dm( iscf )
      use precision
      use units, only: eV
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: qa, lasto, iphorb, iaorb, no_u, no_s, indxuo, 
     &                    qtot, no_l
      use sys, only: die, bye
      use Kpoint_grid
      use m_energies,   only: Ebs, DEHarr, Ecorrec, Ef, Entrop
      use m_rmaxh
      use m_eo
      use m_gamma
      use parallel, only: IONode
      use parallel, only: worker
      use m_pexsi_solver,  only: pexsi_solver
#ifdef MPI
      use mpi_siesta
      use m_mpi_utils,     only:  globalize_sum
#endif
#ifdef CDF
      use iodmhs_netcdf, only: write_dmh_netcdf
#endif
      use m_dminim,        only : dminim
      use m_zminim,        only : zminim
      use m_ordern,        only : ordern
      use m_steps,         only : istp
      use m_spin,          only : nspin, qs, efs
#ifdef TRANSIESTA
      use m_ts_options, only : mixH
      use m_ts_global_vars, only : TSiscf, TSinit, TSrun, Hold
      use m_ts_options,     only : TSmode
      use m_ts_in_siesta,   only : transiesta
#endif /* TRANSIESTA */
      implicit none
!     Input variables
      integer, intent(in) :: iscf
!     Local variables
      logical             :: CallDiagon
      integer             :: io, ispin, nnz
      real(dp)            :: const ! Temporary factor used in
                             ! calculation of Harris forces. e1>e2 to
                             ! signal that we do not want DOS weights
      real(dp), parameter :: e1 = 1.0_dp, e2 = -1.0_dp
#ifdef MPI
      real(dp)            :: buffer1, freeEnergyCorrection
      integer             :: mpierr
#endif

!-------------------------------------------------------------------- BEGIN

      if (worker) call timer( 'compute_dm', 1 )

      call MPI_Bcast(isolve,1,MPI_integer,0,true_MPI_Comm_World,mpierr)

      if (isolve .eq. SOLVE_PEXSI) then
         ! Prototype interface
         call pexsi_solver(no_u, no_l, nspin,
     $              maxnh, numh, listhptr, listh,
     $              H, S, qtot, Dscf, Escf,
     $              ef, freeEnergyCorrection, temp)
         if (ionode) write (6,"(/a,f14.6)")
     &    'Free energy correction term:', freeEnergyCorrection
         Entrop = -freeEnergyCorrection / temp

      endif
      if (.not. worker) RETURN

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
     &              eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
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
        Ecorrec = 0.0_dp
        Entrop = 0.0_dp
      elseif (isolve .eq. SOLVE_MINIM) then
        if (nspin==4) call die('ERROR: Non-collinear spin calculations
     &                          not yet implemented with the OMM!')
        if (gamma) then
          call dminim(.false., PreviousCallDiagon, iscf, istp, no_l,
     &                 nspin, no_u, maxnh, numh, listhptr, listh, Dscf,
     &                 eta, qs, h, s, h_kin)
        else
          call zminim(.false., PreviousCallDiagon, iscf, istp, no_l,
     &                 nspin, no_u, maxnh, numh, listhptr, listh, Dscf,
     &                 eta, qs, no_s, xijo, indxuo, nkpnt, kpoint,
     &                 kweight, h, s, h_kin)
        end if
        Ecorrec = 0.0_dp
        Entrop = 0.0_dp
        PreviousCallDiagon=.false.
#ifdef TRANSIESTA
! TSS Begin
      elseif (TSmode .and. TSinit) then
        call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh,
     &              H, S, qtot, fixspin, qs, temp, e1, e2,
     &              gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     &              occtol, iscf, neigwanted)

        Ecorrec = 0._dp

      else if (TSrun) then
        call timer('TS',1)
        TSiscf=TSiscf+1

        call transiesta(ucell,gamma,lasto,iaorb,
     &       no_l, no_s, no_u, nspin, xa, na_u, maxnh, numh, listhptr,
     &       listh, xijo, indxuo, H, S, Dscf, Escf, Ef,
     &       TSiscf,qtot)
        call timer('TS',2)

        Ecorrec = 0._dp
! TSS End
#endif /* TRANSIESTA */
      else
        !call die('siesta: ERROR: wrong solution method')
      endif

#ifdef TRANSIESTA
! Save present H matrix
          if(mixH) then
             do ispin = 1,nspin
                do io = 1,maxnh
                   Hold(io,ispin) = H(io,ispin)
                enddo
             enddo
          end if
#endif /* TRANSIESTA */

!     Harris-functional energy 
      DEharr = 0.0_dp
      Ebs = 0.0_dp
      do ispin = 1,nspin
!       const factor takes into account that there are two nondiagonal
!       elements in non-collinear spin density matrix, stored as
!       ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
!       ispin=4 => Imag(D12)
        const = 1._dp
        if (ispin .gt. 2) const = 2._dp
        do io = 1,maxnh
          DEharr = DEharr + H(io,ispin) * const * 
     &                 ( Dscf(io,ispin) - Dold(io,ispin) )
          Ebs = Ebs + H(io,ispin) * const * Dscf(io,ispin)
!         Alternative form
!         Ebs = Ebs + S(io) * const * Escf(io,ispin)
        enddo
      enddo
#ifdef MPI
!     Global reduction of DEharr
      call globalize_sum( DEharr, buffer1 )
      DEharr = buffer1
      call globalize_sum(Ebs,buffer1)
      Ebs = buffer1
#endif
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

          if (IOnode .and. iscf==1) then
            call write_orb_indx( na_u, na_s, no_u, no_s, isa, xa,
     .                           iaorb, iphorb, indxuo, nsc, ucell )
          endif

      call timer( 'compute_dm', 2 )

!----------------------------------------------------------------------- END
      END subroutine compute_dm
      END MODULE m_compute_dm
