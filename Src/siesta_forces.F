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
      use parallel,     only : IOnode, SIESTA_worker
      use m_state_analysis
      use m_steps
      use sys, only : die, bye
      use sparse_matrices, only: H, Hold
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance
      use siesta_geom,   only: na_u           ! Number of atoms in unit cell
      use m_energies,    only: Etot           ! Total energy
      use m_forces,      only: fa, cfa        ! Forces and constrained forces
      use m_stress,      only: cstress        ! Constrained stress tensor
      use siesta_master, only: forcesToMaster ! Send forces to master prog
      use siesta_master, only: siesta_server  ! Is siesta a server?
      use write_subs
      use units,         only: eV, Ang
      use m_pexsi_solver,        only: prevDmax
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

      if (SIESTA_worker)  then
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
      endif ! SIESTA_worker

      call MPI_Bcast(last,1,MPI_logical,0,true_MPI_Comm_World,mpierr)
      call MPI_Bcast(iscf,1,MPI_integer,0,true_MPI_Comm_World,mpierr)

      DO WHILE (.not. last)

       if (SIESTA_worker) then
          call timer( 'IterSCF', 1 )
          if (cml_p)
     &         call cmlStartStep( xf=mainXML, type='SCF', index=iscf )
       endif

       if (mixH) then
          if ( first .and. SIESTA_worker) call setup_hamiltonian( iscf )
          call compute_dm( iscf )
          if (SIESTA_worker) call setup_hamiltonian( iscf )
       else 
          if (SIESTA_worker) call setup_hamiltonian( iscf )
          call compute_dm( iscf )
       end if                   ! mixH
       if (SIESTA_worker) call mixer( iscf, dDmax )

       if (SIESTA_worker) then
          prevDmax = dDMax
          call scfconvergence_test( first, iscf, last, dDmax,
     &                            conv_harris, conv_etot )
       endif ! SIESTA_worker

       call MPI_Bcast(last,1,MPI_logical,0,true_MPI_Comm_World,mpierr)

       if (SIESTA_worker) then
!       At this point last is true if SCF has converged. But this
!       changed below
        SCFconverged = last

        if (monitor_forces_in_scf) then
           allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
           fa_old = fa
           if (mixH) Hold(:,:) = H(:,:)
           call final_H_f_stress(iscf)
           if (mixH) H(:,:) = Hold(:,:)     
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
       endif ! SIESTA_worker

        iscf = iscf + 1
        first = .false.

        if (SIESTA_worker) then
!       The last scheduled iteration is used to compute forces
!       and stresses, hence the exit if iscf==nscf.
!       The ">" clause catches the nscf==1 case. For this,
        if (iscf < min_nscf ) then
           last = .false.
        else if (iscf >= nscf) then
           last = .true.
        endif
        endif ! SIESTA_worker
        call MPI_Bcast(last,1,MPI_logical,0,true_MPI_Comm_World,mpierr)
      ENDDO

      if (.not. SIESTA_worker) RETURN

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

      ! If siesta is running as a subroutine, send forces to master program
      if (siesta_server) 
     &  call forcesToMaster( na_u, Etot, cfa, cstress )


!------------------------------------------------------------------------ END
      END subroutine siesta_forces
      END module m_siesta_forces
