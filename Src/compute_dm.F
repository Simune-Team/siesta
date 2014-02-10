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
      use m_energies,   only: Ebs, Ecorrec, Ef, Entropy
      use m_rmaxh
      use m_eo
      use m_gamma
      use parallel, only: IONode
      use parallel, only: SIESTA_worker
      use m_compute_ebs_shift, only: compute_ebs_shift
      use m_pexsi_solver,  only: pexsi_solver
      use m_hsx, only      : write_hs_formatted
c$$$      use m_matrixwrite
      use mpi_siesta
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
      use m_ts_global_vars, only : TSiscf, TSinit, TSrun
      use m_ts_options,     only : TSmode
      use m_ts_in_siesta,   only : transiesta
#endif /* TRANSIESTA */
      implicit none
!     Input variables
      integer, intent(in) :: iscf

      real(dp)            :: delta_Ebs, delta_Ef
      logical             :: CallDiagon
      integer             :: nnz

      ! e1>e2 to signal that we do not want DOS weights
      real(dp), parameter :: e1 = 1.0_dp, e2 = -1.0_dp
      real(dp)            :: buffer1
      integer             :: mpierr

!       character(15)            :: filename, indexstr
!       character(15), parameter :: fnameform = '(A,A,A)'


!-------------------------------------------------------------------- BEGIN

      if (SIESTA_worker) call timer( 'compute_dm', 1 )

      call MPI_Bcast(isolve,1,MPI_integer,0,true_MPI_Comm_World,mpierr)

      if (SIESTA_worker) then
         ! Save present density matrix 
         Dold(1:maxnh,1:nspin) = Dscf(1:maxnh,1:nspin)
         Eold(1:maxnh,1:nspin) = Escf(1:maxnh,1:nspin)
      endif

      ! Compute shift in Tr(H*DM) for fermi-level bracketting
      ! Use the current H, the previous iteration H, and the
      ! previous iteration DM

      if (SIESTA_worker) then
         if (iscf > 1) then
            call compute_Ebs_shift(Dscf,H,Hprev,delta_Ebs)
            delta_Ef = delta_Ebs / qtot
            if (ionode) then
               write(6,"(a,f16.5)")
     $              "Estimated change in band-structure energy:",
     $              delta_Ebs/eV, "Estimated shift in E_fermi: ",
     $              delta_Ef/eV
            endif
               
         else
            delta_Ebs = 0.0_dp
            delta_Ef  = 0.0_dp
         endif
      endif

      if (isolve .eq. SOLVE_PEXSI) then
         ! Prototype interface
         call pexsi_solver(iscf, no_u, no_l, nspin,
     $              maxnh, numh, listhptr, listh,
     $              H, S, qtot, Dscf, Escf,
     $              ef, Entropy, temp, delta_Ef)
         if (ionode) write (6,"(/a,f14.6)") 'Entropy/k:', Entropy
      endif
      if (.not. SIESTA_worker) RETURN

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

      if (isolve .eq. MATRIX_WRITE) then
!             write(indexstr,'(I15)') iscf
!             write(filename,fnameform) 'H_', trim(adjustl(indexstr)), 
!      &                                '.matrix'
!             call write_global_matrix( no_s, no_l, maxnh, numh, listh,
!      &           H(1:maxnh,1), filename )
!  
!             write(filename,fnameform) 'S_', trim(adjustl(indexstr)), 
!      &                                '.matrix'

!        Note: only one-shot for now
         call write_hs_formatted(no_u, nspin,
     $        maxnh, numh, listhptr, listh, H, S)
         call bye("End of run after writing H.matrix and S.matrix")

c$$$        call write_global_matrix_singlenodewrite( 
c$$$     &           no_u, no_s, maxnh, numh, listhptr, listh, 
c$$$     &           H(:,1), 'H.matrix')
c$$$
c$$$        call write_global_matrix_singlenodewrite( 
c$$$     &           no_u, no_s, maxnh, numh, listhptr, listh, 
c$$$     &           S, 'S.matrix')

      elseif ((isolve .eq. SOLVE_DIAGON) .or. (CallDiagon)) then
        call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh, 
     &              H, S, qtot, fixspin, qs, temp, e1, e2,
     &              gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &              occtol, iscf, neigwanted)
        Ecorrec = 0.0_dp
        PreviousCallDiagon=.true.
        if (ionode) write (6,"(/a,f14.6)") 'Entropy/k:', Entropy
      elseif (isolve .eq. SOLVE_ORDERN) then
        if (.not. gamma) call die("Cannot do O(N) with k-points.")
        call ordern(usesavelwf, ioptlwf, na_u, no_u, no_l, lasto,
     &               isa, qa, rcoor, rmaxh, ucell, xa, iscf,
     &               istp, ncgmax, etol, eta, qtot, maxnh, numh,
     &               listhptr, listh, H, S, chebef, noeta, rcoorcp,
     &               beta, pmax, Dscf, Escf, Ecorrec, nspin, qs )
        Ecorrec = 0.0_dp
        Entropy = 0.0_dp
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
        Entropy = 0.0_dp
        PreviousCallDiagon=.false.
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

          if (IOnode .and. iscf==1) then
            call write_orb_indx( na_u, na_s, no_u, no_s, isa, xa,
     .                           iaorb, iphorb, indxuo, nsc, ucell )
          endif

!     Normalize density matrix to exact charge
!     Placed here for now to avoid disturbing EHarris

#ifdef TRANSIESTA      
         if (.not.TSrun ) call normalize_dm( first= .false. )
#else
         call normalize_dm(first=.false.)
#endif  /* TRANSIESTA */

! Save Hamiltonian used to generate the DM in this iteration
      Hprev = H

      call timer( 'compute_dm', 2 )
      if (ionode) call memory_snapshot("after compute_DM")

!----------------------------------------------------------------------- END
      END subroutine compute_dm
      END MODULE m_compute_dm
