      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)

      use precision,             only: dp
      use siesta_cml,            only: cml_p, mainXML
      use siesta_cml,            only: cmlStartStep, cmlEndStep
      use sparse_matrices,       only: H, Hold, Dscf, Dold
      use m_compute_max_diff,    only: compute_max_diff
      use m_state_init,          only: state_init
      use m_setup_hamiltonian,   only: setup_hamiltonian
      use m_setup_H0,            only: setup_H0
      use m_compute_dm,          only: compute_dm
      use m_EKS_Harris,          only: compute_EKS_Harris
      use m_check_scf_convergence, only: check_scf_convergence
      use m_post_scf_work,       only: post_scf_work
      use m_mixer,               only: mixer
      use siesta_options,        only: require_harris_convergence
      use siesta_options,        only: require_energy_convergence
      use siesta_options,        only: energy_tolerance
      use siesta_options,        only: harris_tolerance
      use siesta_options,        only: SCFMustConverge, nscf, mixH
      use siesta_options,        only: mix_after_convergence
      use parallel,              only: IOnode
      use m_state_analysis,      only: state_analysis
      use m_steps,               only: inicoor
      use sys,                   only: die, bye
      use m_convergence,         only: converger_t
      use m_convergence,         only: reset, set_tolerance
      use parallel,              only: IONode
#ifdef TRANSIESTA
      use m_ts_global_vars,      only: ts_istep
#endif /* TRANSIESTA */

      implicit none

      integer, intent(in)  :: istep

      integer            :: iscf
      logical            :: last, SCFconverged
      type(converger_t)  :: conv_harris, conv_etot
      real(dp)           :: E_Harris, E_KS_Good
      real(dp)           :: dDmax, dHmax

#ifdef TRANSIESTA
! Used in force constant runs
      ts_istep = istep
#endif /* TRANSIESTA */
      
!---------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( 'PRE siesta_forces' )
#endif
!     Initialization tasks for a given geometry, including
!     the determination of DM_in
      call state_init( istep )

!     The current structure of the loop tries to reproduce the
!     historical Siesta usage. It should be made more clear.

!     Start of SCF loop
      iscf  = 1
      last  = .false.
      SCFconverged = .false.
      if (require_harris_convergence) then
        call reset(conv_harris)
        call set_tolerance(conv_harris,harris_tolerance)
      endif
      if (require_energy_convergence) then
        call reset(conv_etot)
        call set_tolerance(conv_etot,energy_tolerance)
      endif

!     The SCF loop is executed at least once, 
!     except if nscf <= 0.
      if (nscf <= 0) then
         last = .true.
      endif

      call setup_H0()

      ! Beginning of SCF loop
      ! Starting point is H_in
      call timer( 'IterSCF', 1 )
      if (cml_p)
     &     call cmlStartStep( xf=mainXML, type='SCF', index=iscf )
      ! Note that, if mixing H, H_in could also be read from file ...
      H(:,:) = 0.0_dp
      call setup_hamiltonian( iscf ) ! Gets E_KS(DM_in) (non-variational) and H_in
      call compute_max_diff(Hold,H,dHmax)

      DO WHILE (.not. last)

        ! The underlying reason for using "last" is to allow
        ! TranSiesta to keep looping after the initial diagon-style
        ! initialization has converged.
        ! This should be made more clear.


        call compute_dm( iscf ) ! Gets DM_out (Dscf) and keeps DM_in (Dold)
        call compute_max_diff(Dold,Dscf,dDmax)
        call compute_EKS_Harris( E_Harris, E_KS_Good )

        ! Checks convergence

        call check_scf_convergence( iscf, last, E_Harris, E_KS_Good,
     &                            dDmax, dHmax, conv_harris, conv_etot )
        SCFconverged = last
        ! In TranSiesta runs, last is re-set to .false. to switch
        ! to GF mode (to be cleaned up)

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. (iscf==1)) then
           call timer( 'IterSCF', 3 )
        endif
        if (cml_p) call cmlEndStep(mainXML)


        if (.not. SCFconverged) then
           call timer( 'IterSCF', 1 )
           if (cml_p) call cmlStartStep(mainXML)
           if (mixH) then
              call setup_hamiltonian( iscf )  ! Gets H_out (H) and keeps H_in (Hold)
              call mixer (iscf)               ! Mixes H_in and H_out to give H_in_next
           else
              call mixer (iscf)               ! Mixes DM_in and DM_out to give DM_in_next
              call setup_hamiltonian (iscf)   ! Gets H_in_next from DM_in_next, and E_KS(DM_in)
           endif
           call compute_max_diff(Hold,H,dHmax)
        endif

        ! Save DM and H here for restarts

        iscf = iscf + 1

!       The last scheduled iteration is used to compute forces
!       and stresses, hence the exit if iscf==nscf.
!       The ">" clause catches the nscf==1 case. For this,
        if (iscf >= nscf) then
          last = .true.
        endif

      ENDDO

      !     End of standard SCF loop.

      if ((SCFMustConverge).and.(.not.SCFconverged)) then
        call die('SCF did not converge in maximum number of steps.')
      endif

      ! Do one more pass to compute forces and stresses (and energies)
      ! with DM_out coming out of the cycle

      call post_scf_work( iscf )

      ! Perform any other analyses of the electronic structure
      ! for this geometry
      call state_analysis( istep )

!--------------------------------------------------------------------------- END
      END subroutine siesta_forces
      END module m_siesta_forces
