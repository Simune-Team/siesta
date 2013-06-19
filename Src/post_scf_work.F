      MODULE m_post_scf_work
      private
      public :: post_scf_work

      CONTAINS

      subroutine post_scf_work( iscf )
      USE m_final_H_f_stress, only: final_H_f_stress
      USE siesta_options
      use sparse_matrices
      use siesta_geom
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

      integer  :: iscf
      external :: plcharge

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE post_scf_work' )
#endif
      call timer( 'PostSCF', 1 )
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
     &                eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
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

!     If converged, make one last iteration to find forces and
!     stress
!     Note that here we are using the Dscf guessed by mixer,
!     not the last Dscf really computed from Wavefunctions
      call final_H_f_stress( iscf )

!     Print energies
      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec +
     &       Emad + Emm + Emeta
      Eharrs = Etot + DEharr
      FreeE  = Etot - Temp * Entropy

!
!     WHY IS THIS?? To make Dscf (which might be re-used in the next
!     geometry iteration) equal to the last DM used to construct the
!     hamiltonian, instead of the output DM.
!
!       We need the last D really computed from Wavefunctions,
!       not the one guessed by mixer as the D_in for a hypothetical
!       next iteration. BUT in the old version of the program
!       Dold was set equal to Dscf at the end of setup_hamiltonian,
!       also in the "last" iteration.  
      
      Dscf(1:maxnh,1:nspin) = Dold(1:maxnh,1:nspin)
      Escf(1:maxnh,1:nspin) = Eold(1:maxnh,1:nspin)

      if (dumpcharge) then
        call plcharge( no_s, na_s, no_u, maxnh, maxna, nspin,
     &       isa, iphorb, indxuo, lasto,
     &       scell, nsc, xa, rmaxo, datm )
      endif

      call timer( 'PostSCF', 2 )
#ifdef DEBUG
      call write_debug( '  POS post_scf_work' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine post_scf_work

      END MODULE m_post_scf_work
