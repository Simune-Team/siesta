      MODULE m_post_scf_work
      private
      public :: post_scf_work

      CONTAINS

      subroutine post_scf_work( iscf )

      USE m_final_H_f_stress, only: final_H_f_stress
      USE siesta_options
      use sparse_matrices, only: Dscf, Escf, Dold, Eold, maxnh
      use siesta_geom
      use atomlist, only: lasto, rmaxo, datm, indxuo, no_s, no_u, iphorb
      use m_energies
      use neighbour,   only: maxna=>maxnna   ! For plcharge...
      use m_spin,         only: nspin
      implicit none

      integer, intent(in)  :: iscf

      integer  :: io, is
      external :: plcharge

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE post_scf_work' )
#endif
      call timer( 'PostSCF', 1 )
!     If converged, make one last iteration to find forces and
!     stress
!     Note that, if mixing the DM, we are using here
!     the Dscf guessed by mixer for a hypothetical next iteration,
!     not the last Dscf really computed from Wavefunctions
!     If mixing the Hamiltonian, we are using the output DM from
!     'compute_dm'.
!     We are using the Escf coming out of compute_dm also.

      ! This call might also involve "grid-cell-sampling" to
      ! refine the grid contribution to the energy and forces/stress.
      call final_H_f_stress( iscf )

!     Print energies
      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec +
     &       Emad + Emm + Emeta
      ! NOTE: DEharr is only computed by compute_dm, and 
      ! corresponds to Tr(H*(DM_out-DM_in))
      ! Since this Etot is computed with DM_next, this value
      ! for Eharrs is suspect. To make things worse, this value
      ! is never actually printed, since printing routines use
      ! Eharrs1, which is the last value of Eharrs computed
      ! in the scf cycle (set in scfconvergence_test).
      ! Apparently, Eharrs1 is saved in case a further grid-cell
      ! sampling is carried out, presumably to avoid recomputing
      ! it with an "impure" H. In fact, Eharrs1 might be what we
      ! need everywhere, particularly in the computation of the
      ! FreeEHarris value used if performing Harris-oriented runs,
      ! as in basis-optimization (the only magnitude depending 
      ! explicitly on Eharrs, and thus exhibiting a mismatch with
      ! what would be obtained from Eharrs1).
      ! From this point of view, what we need is to get rid of this
      ! statement and replace Eharris1 with Eharrs elsewhere.
      ! Note also in general that the "Harris forces" are computed
      ! by the above call to final_H_f_stress, which is using the
      ! wrong DM! This would explain the non-conservation of the energy
      ! in some Harris-based MD runs.
      ! All this gives further support to the idea that the DM coming
      ! out of the SCF cycle for post-processing should be the last
      ! DM_out...

      ! Eharrs = Etot + DEharr

      !
      FreeE  = Etot - Temp * Entropy
!
! The following block has been moved here from 
! setup_H_f_stress for clarity

! If mixing H, Dold is the last DM used as input in setup_hamiltonian,
! and Dscf is the output DM from compute_dm. So they are different.
! The following block copies them.

! If mixing the DM, Dold and Dscf are equal, and have been set by mixer
! to the "predicted input DM" for a hypothetical next iteration.
! So the following block is superfluous.

! Save present density matrix ........................................
      do is = 1,nspin
        do io = 1,maxnh
          Dold(io,is) = Dscf(io,is)
          Eold(io,is) = Escf(io,is)
        enddo
      enddo        

!   NOTE that this block will cancel the effect of the previous one...
!   but before removing them we should understand this more deeply.
!
!     WHY IS THIS?? To make Dscf (which might be re-used in the next
!     geometry iteration) equal to the last DM used to construct the
!     hamiltonian, instead of the output DM.
!     (in the absence of the above block...)
!
!       We need (DO WE?) the last D really computed from Wavefunctions,
!       not the one guessed by mixer as the D_in for a hypothetical
!       next iteration. BUT in the old version of the program
!       Dold was set equal to Dscf at the end of setup_hamiltonian,
!       also in the "last" iteration.  
      
      Dscf(1:maxnh,1:nspin) = Dold(1:maxnh,1:nspin)
      Escf(1:maxnh,1:nspin) = Eold(1:maxnh,1:nspin)

!
!
!
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
