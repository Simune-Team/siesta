      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      use siesta_cml
      use m_state_init
      use m_setup_hamiltonian
      use m_final_H_f_stress,    only: final_H_f_stress
      use m_setup_H0
      use m_compute_dm
      use m_compute_max_diff
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use parallel,     only : IOnode
      use m_state_analysis
      use m_steps
      use sys, only : die, bye
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
      use sparse_matrices, only: H, Hold, Dold, Dscf
      use m_save_density_matrix, only: save_density_matrix
#ifdef TRANSIESTA
      use m_ts_options, only : mixH
      use m_ts_global_vars, only : ts_istep
#endif /* TRANSIESTA */
      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: first, last, SCFconverged
      real(dp) :: dmax      ! Max. change in DM or H elements
      type(converger_t)  ::  conv_harris, conv_etot
      real(dp), allocatable  :: fa_old(:,:), Hsave(:,:)
#ifdef TRANSIESTA
! Used in force constant runs
      ts_istep = istep
#endif /* TRANSIESTA */
      
!---------------------------------------------------------------------- BEGIN

!     Initialization tasks for a given geometry
      call state_init( istep )

!     The current structure of the loop tries to reproduce the
!     historical Siesta usage. It should be made more clear.

!     Start of SCF loop
      iscf  = 1
      first = .true.
      last  = .false.

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
        call timer( 'IterSCF', 1 )
        if (cml_p)
     &         call cmlStartStep( xf=mainXML, type='SCF', index=iscf )
       
        if (mixH) then
           if ( first) call setup_hamiltonian( iscf )
           call compute_dm( iscf )
           call setup_hamiltonian( iscf )
           call compute_max_diff(Hold,H,dmax)
        else 
           call setup_hamiltonian( iscf )
           call compute_dm( iscf )
           call compute_max_diff(Dold,Dscf,dmax)
        end if 

        call scfconvergence_test( first, iscf, last, dmax,
     &                            conv_harris, conv_etot,
     &                            SCFconverged )

        if (monitor_forces_in_scf) then
           allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
           allocate(Hsave(size(H,dim=1),size(H,dim=2)))
           fa_old = fa
           Hsave(:,:) = H(:,:)
           call final_H_f_stress(iscf)
           H(:,:) = Hsave(:,:)
           if (ionode) then
              print *, "Max diff in force (eV/Ang): ",
     $                 maxval(abs(fa-fa_old))*Ang/eV
           endif
           deallocate(fa_old,Hsave)
           call siesta_write_forces()
           call siesta_write_stress_pressure()
        endif

        if ((.not. SCFconverged) .OR. mix_after_convergence) then
           ! Mix for next step
           call mixer( iscf )
           ! Save for possible restarts
           ! Could save H in the future
           call save_density_matrix("mixed")
        endif

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )

        iscf = iscf + 1
        first = .false.

!       The last scheduled iteration is used to compute forces
!       and stresses, hence the exit if iscf==nscf.
!       The ">" clause catches the nscf==1 case. For this,
        if (iscf < min_nscf ) then
           last = .false.
        else if (iscf >= nscf) then
           last = .true.
        endif
        if (cml_p) call cmlEndStep(mainXML)
      ENDDO

!     If we have not converged, but get to here, we are out
!     of steps. We may want to error out.
      if ((SCFMustConverge).and.(.not.SCFconverged)) then
        call die('SCF did not converge in maximum number of steps.')
      endif

      ! Save final DM_out for further analysis offline
      if (.not. mix_after_convergence) then
          call save_density_matrix("dmout")
          if (MixH) then
             ! We could make sure that we keep the H actually used 
             ! to generate the last DM, if needed.
             !!!    H = Hold
          endif
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
