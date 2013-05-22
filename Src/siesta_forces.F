      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      use mpi_siesta
      use siesta_cml
      use m_state_init
      use m_setup_hamiltonian
      use m_final_H_f_stress,    only: final_H_f_stress
      use m_setup_H0
      use m_compute_dm
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use parallel,     only : IOnode, worker
      use m_state_analysis
      use m_steps
      use sys, only : die, bye
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance
      use write_subs
      use m_forces,              only: fa
      use units,                 only: eV, Ang
#ifdef TRANSIESTA
      use m_ts_options, only : mixH
      use m_ts_global_vars, only : ts_istep
#endif /* TRANSIESTA */
      implicit none
      integer, intent(in)  :: istep
      integer  :: iscf
      logical  :: first, last, SCFconverged
      real(dp) :: dDmax     ! Max. change in density matrix elements
      type(converger_t)  ::  conv_harris, conv_etot
      real(dp), allocatable  :: fa_old(:,:)
#ifdef TRANSIESTA
! Used in force constant runs
      ts_istep = istep
#endif /* TRANSIESTA */

      integer mpierr
!---------------------------------------------------------------------- BEGIN

      if (worker)  then
!     Initialization tasks for a given geometry
      call state_init( istep )

!     The current structure of the loop tries to reproduce the
!     historical Siesta usage. It should be made more clear.

!     Start of SCF loop
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

!     The SCF loop is executed at least once, 
!     except if nscf <= 0.
      if (nscf <= 0) then
         last = .true.
      endif

      call setup_H0()
      endif ! worker

      call MPI_Bcast(last,1,MPI_logical,0,true_MPI_Comm_World,mpierr)
      call MPI_Bcast(iscf,1,MPI_integer,0,true_MPI_Comm_World,mpierr)

      DO WHILE (.not. last)

       if (worker)  call timer( 'IterSCF', 1 )

#ifndef TRANSIESTA
       if (worker) then
        if (cml_p) 
     &    call cmlStartStep(xf=mainXML, type='SCF', index=iscf)

!       Will use H0 computed above
        call setup_hamiltonian( iscf )
!!        call printmemory( 6, 0 )

        if (h_setup_only) then
           call timer( 'all', 2 )  ! New call to close the tree
           call timer( 'all', 3 )
           call bye("H-Setup-Only requested")
           STOP
        endif
       endif  ! worker

       call compute_dm( iscf )

       if (worker) call mixer( iscf, dDmax )

#else /* TRANSIESTA */
       if (worker) then
        if (cml_p)
     &      call cmlStartStep( xf=mainXML, type='SCF', index=iscf )

! Adding the possibility of Mixing H.

        if (mixH) then
          if ( first ) call setup_hamiltonian( iscf )
          call compute_dm( iscf )
          call setup_hamiltonian( iscf )
        else 
          call setup_hamiltonian( iscf )
          call compute_dm( iscf )
        end if ! mixH
      	call mixer( iscf, dDmax )
       endif ! worker
#endif /* TRANSIESTA */

       if (worker) then
          call scfconvergence_test( first, iscf, last, dDmax,
     &                            conv_harris, conv_etot )
       endif ! worker

       call MPI_Bcast(last,1,MPI_logical,0,true_MPI_Comm_World,mpierr)

       if (worker) then
!       At this point last is true if SCF has converged. But this
!       changed below
        SCFconverged = last

        if (monitor_forces_in_scf) then
           allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
           fa_old = fa
           call final_H_f_stress(iscf)
           if (ionode) then
              print *, "Max diff in force (eV/Ang): ",
     $                 maxval(abs(fa-fa_old))*Ang/eV
           endif
           deallocate(fa_old)
           call siesta_write_forces()
           call siesta_write_stress_pressure()
        endif

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )
        if (cml_p) call cmlEndStep(mainXML)
       endif ! worker

        iscf = iscf + 1
        first = .false.

        if (worker) then
!       The last scheduled iteration is used to compute forces
!       and stresses, hence the exit if iscf==nscf.
!       The ">" clause catches the nscf==1 case. For this,
        if (iscf < min_nscf ) then
           last = .false.
        else if (iscf >= nscf) then
           last = .true.
        endif
        endif ! worker
        call MPI_Bcast(last,1,MPI_logical,0,true_MPI_Comm_World,mpierr)
      ENDDO

      if (.not. worker) RETURN

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
