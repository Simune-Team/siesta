      module m_siesta_forces
      private
      
      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      use siesta_cml
      use m_state_init
      use m_setup_hamiltonian
      use m_setup_H0
      use m_compute_dm
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use parallel,     only : IOnode
      use m_state_analysis
      use m_steps
      use sys, only : die
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance

      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: first, last, SCFconverged
      real(dp) :: dDmax     ! Max. change in density matrix elements
      type(converger_t)  ::  conv_harris, conv_etot
      
!---------------------------------------------------------------------- BEGIN
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 10 )
#endif
C     Initialization tasks for a given geometry
#ifdef DEBUG
      call write_debug( '  PRE state_init' )
#endif
      call state_init( istep )

C     The current structure of the loop tries to reproduce the
C     historical Siesta usage. It should be made more clear.

C     Start of SCF loop
      iscf  = 1
      first = .true.
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

C     The SCF loop is executed at least once, 
C     except if nscf <= 0.
      if (nscf <= 0) then
         last = .true.
      endif

#ifdef DEBUG
        call write_debug( '  PRE setup_H0' )
#endif
      call setup_H0()

      DO WHILE (.not. last)
        call timer( 'IterSCF', 1 )
#ifdef DEBUG
        call write_debug( '  PRE setup_hamiltonian' )
#endif

          if (cml_p) 
     $         call cmlStartStep(xf=mainXML, type='SCF', index=iscf)
C       Will use H0 computed above
        call setup_hamiltonian( iscf )

#ifdef DEBUG
        call write_debug( '  PRE compute_dm' )
#endif
        call timer( 'compute_dm', 1 )
        call compute_dm( iscf )
        call timer( 'compute_dm', 2 )

#ifdef DEBUG
        call write_debug( '  PRE mixer' )
#endif

        call mixer( iscf, dDmax )

#ifdef DEBUG
        call write_debug( '  PRE scfconvergence_test' )
#endif
        call scfconvergence_test( first, iscf, last, dDmax,
     $                              conv_harris, conv_etot )
        ! At this point last is true if SCF has converged. But this
        ! changed below
        SCFconverged = last

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )

        iscf = iscf + 1
        first = .false.

C       The last scheduled iteration is used to compute forces
C       and stresses, hence the exit if iscf==nscf.
C       The ">" clause catches the nscf==1 case. For this,
        if (iscf >= nscf) then
          last = .true.
        endif
        if (cml_p) call cmlEndStep(mainXML)
      ENDDO

        ! If we have not converged, but get to here, we are out
        ! of steps. We may want to error out.
        if ((SCFMustConverge).and.(.not.SCFconverged)) then
          call die('SCF did not converge in maximum number of steps.')
        endif

C       End of standard SCF loop.
C       Do one more pass to compute forces and stresses
      if (last) then
#ifdef DEBUG
        call write_debug( '  PRE post_scf_work' )
#endif
        call post_scf_work( iscf )
      else
        STOP 'Bad logic in SCF loop'
      endif

#ifdef DEBUG
        call write_debug( '  PRE state_analysis' )
#endif
      call state_analysis( istep )

#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 0 )
#endif
!--------------------------------------------------------------------------- END
        END subroutine siesta_forces
      END module m_siesta_forces
