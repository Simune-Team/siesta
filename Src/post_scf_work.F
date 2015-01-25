      MODULE m_post_scf_work
      private
      public :: post_scf_work

      CONTAINS

      subroutine post_scf_work( iscf )

      USE m_final_H_f_stress, only: final_H_f_stress
      USE siesta_options
      use sparse_matrices, only: Dscf, Escf, Dold, maxnh
      use sparse_matrices, only: H, S
      use sparse_matrices, only: xijo, listh, listhptr, numh
      use siesta_geom
      use parallel, only: ionode
      use atomlist, only: lasto, rmaxo, datm, indxuo, no_s, no_u,
     &                    iphorb, no_l, qtot
      use m_energies
      use neighbour,   only: maxna=>maxnna   ! For plcharge...
      use m_spin,         only: nspin, qs, efs
      use m_dminim,       only : dminim
      use m_zminim,       only : zminim
      use m_steps,        only : istp
      use m_compute_dm,   only : PreviousCallDiagon
      use m_eo
      use Kpoint_grid
      use m_gamma
      implicit none

      integer, intent(in)  :: iscf

      integer  :: io, is
      external :: plcharge

      call timer( 'PostSCF', 1 )
!
!     If converged, make one last iteration to find forces and stress

!     If we use the minimization routine, the energy-density
!     matrix is not calculated inside the SCF step, and
!     so this must be done now

      if ((isolve .eq. SOLVE_MINIM) .and.
     &    (.not. PreviousCallDiagon)) then
        if (minim_calc_eigenvalues) then
          call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     &                numh, listhptr, listh, numh, listhptr, listh, 
     &                H, S, qtot, fixspin, qs, temp, 1.0_dp, -1.0_dp,
     &                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &                eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &                occtol, iscf, neigwanted)
          Ecorrec = 0.0_dp
        else
          if (gamma) then
            call dminim(.true., .false., iscf, istp, no_l, nspin, no_u,
     &                  maxnh, numh, listhptr, listh, Escf, eta, qs)
          else
            call zminim(.true., .false., iscf, istp, no_l, nspin, no_u,
     &                  maxnh, numh, listhptr, listh, Escf, eta, qs,
     &                  no_s, xijo, indxuo, nkpnt, kpoint, kweight)
          end if
        endif
      endif

!     Find final energy, forces and stress.

!     Note that here we are using the Dscf coming out of the
!     scf cycle, which could be DM_out or DM_mixed, depending
!     on whether we are mixing H or DM, and on whether we
!     mix one final time after convergence

      if (ionode) then
         if (mix_after_convergence .and. (.not. MixH)
     $                             .and. (.not. mix_charge)) then
            write(6,"(/,a)")
     $      "Using DM_mixed to compute the final energy and forces"
            write(6,"(a)")
     $      "Consider using 'SCF.MixAfterConvergence F', mixing H" //
     $           " or mixing the charge density"
         else
            write(6,"(/,a)")
     $      "Using DM_out to compute the final energy and forces"
         endif
      endif

      call final_H_f_stress( iscf )

      ! These energies are computed with a different DM from those in the
      ! last SCF step, except if mixing the Hamiltonian, in which case they
      ! are the same.

      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec +
     &       Emad + Emm + Emeta

!     Since this Etot is computed with DM_mixed (the DM predicted by
!     'mixer' for a hypothetical next SCF iteration) or with DM_out
!     (depending on whether a final mixing is done at the end of the
!     cycle) this value for Eharrs is suspect. This value was never
!     actually printed, since printing routines used 'Eharrs1', which is
!     the last value of Eharrs computed in the SCF cycle (set in
!     'scfconvergence_test'). To reproduce that behavior, the following
!     statement has been removed, and Eharris1 replaced by Eharrs
!     everywhere.

      !! Eharrs = Etot + DEharr
      !
      ! Update the Free electronic energy with the new Etot
      ! and the latest entropy. This is not completely 
      ! consistent if using a mixed DM to compute Etot, but
      ! is a reasonable approximation
      !
      FreeE  = Etot - Temp * Entropy
!
      if (dumpcharge) then
        call plcharge( no_s, na_s, no_u, maxnh, maxna, nspin,
     &       isa, iphorb, indxuo, lasto,
     &       scell, nsc, xa, rmaxo, datm )
      endif

      call timer( 'PostSCF', 2 )

      END subroutine post_scf_work

      END MODULE m_post_scf_work
