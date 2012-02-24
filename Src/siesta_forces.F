      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)

      use precision,             only: dp
      use siesta_cml,            only: cml_p, mainXML
      use siesta_cml,            only: cmlStartStep, cmlEndStep
      use m_state_init,          only: state_init
      use m_setup_hamiltonian,   only: setup_hamiltonian
      use m_setup_H0,            only: setup_H0
      use m_compute_dm,          only: compute_dm
      use m_scfconvergence_test, only: scfconvergence_test
      use m_post_scf_work,       only: post_scf_work
      use m_mixer,               only: mixer
      use siesta_options,        only: require_harris_convergence
      use siesta_options,        only: require_energy_convergence
      use siesta_options,        only: energy_tolerance
      use siesta_options,        only: harris_tolerance
      use siesta_options,        only: SCFMustConverge, nscf, mixH
      use parallel,              only: IOnode
      use m_state_analysis,      only: state_analysis
      use m_steps,               only: inicoor
      use sys,                   only: die, bye
      use m_convergence,         only: converger_t
      use m_convergence,         only: reset, set_tolerance
#ifdef TRANSIESTA
      use m_ts_global_vars,      only: ts_istep
#endif /* TRANSIESTA */

      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: last, SCFconverged
      type(converger_t)  ::  conv_harris, conv_etot

#ifdef TRANSIESTA
! Used in force constant runs
      ts_istep = istep
#endif /* TRANSIESTA */
      
!---------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( 'PRE siesta_forces' )
#endif
!     Initialization tasks for a given geometry
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

      DO WHILE (.not. last)

        ! The underlying reason for using "last" is to allow
        ! TranSiesta to keep looping after the initial diagon-style
        ! initialization has converged.
        ! This should be made more clear.

        call timer( 'IterSCF', 1 )
        if (cml_p)
     &      call cmlStartStep( xf=mainXML, type='SCF', index=iscf )

        if (mixH) then
          if ( iscf == 1) call setup_hamiltonian( iscf )
          call compute_dm( iscf )
          call setup_hamiltonian( iscf )
        else 
          call setup_hamiltonian( iscf )
          call compute_dm( iscf )
        end if ! mixH

        call scfconvergence_test( iscf, last,
     &                            conv_harris, conv_etot )
!       At this point 'last' is true if SCF has converged. But this
!       changed below
        SCFconverged = last

      	if (.not. SCFconverged) call mixer( iscf )

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. (iscf==1)) then
           call timer( 'IterSCF', 3 )
        endif

        iscf = iscf + 1

!       The last scheduled iteration is used to compute forces
!       and stresses, hence the exit if iscf==nscf.
!       The ">" clause catches the nscf==1 case. For this,
        if (iscf >= nscf) then
          last = .true.
        endif
        if (cml_p) call cmlEndStep(mainXML)
      ENDDO

!     If we have not converged, but get to here, we are out
!     of steps. We may want to error out.
      if ((SCFMustConverge).and.(.not.SCFconverged)) then
        call die('SCF did not converge in maximum number of steps.')
      endif

!     End of standard SCF loop.
!     Do one more pass to compute forces and stresses
      if (last) then
        call post_scf_work( iscf )
      else
        STOP 'Bad logic in SCF loop'
      endif

      call state_analysis( istep )
#ifdef DEBUG
      call write_debug( 'POS siesta_forces' )
#endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_forces
      END module m_siesta_forces
