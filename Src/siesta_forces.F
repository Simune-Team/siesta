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
      use m_final_H_f_stress,    only: final_H_f_stress
      use m_setup_H0,            only: setup_H0
      use m_compute_dm,          only: compute_dm
      use m_EKS_Harris,          only: compute_EKS_Harris
      use m_check_scf_convergence, only: check_scf_convergence
      use m_post_scf_work,       only: post_scf_work
      use m_mixer,               only: mixer
      use siesta_options,        only: require_harris_convergence
      use siesta_options,        only: require_energy_convergence
      use siesta_options,        only: require_free_energy_convergence
      use siesta_options,        only: energy_tolerance
      use siesta_options,        only: harris_tolerance
      use siesta_options,        only: freeEnergy_tolerance
      use siesta_options,        only: SCFMustConverge, nscf, mixH
      use siesta_options,        only: mix_after_convergence
      use siesta_options,        only: monitor_forces_in_scf
      use m_state_analysis,      only: state_analysis
      use m_steps,               only: inicoor
      use sys,                   only: die, bye
      use m_convergence,         only: converger_t
      use m_convergence,         only: reset, set_tolerance
      use parallel,              only: IONode
      use m_save_density_matrix, only: save_density_matrix
      use m_forces,              only: fa
      use units,                 only: eV, Ang
      use write_subs,            only: siesta_write_forces
      use write_subs,            only: siesta_write_stress_pressure
#ifdef TRANSIESTA
      use m_ts_global_vars,      only: ts_istep
#endif /* TRANSIESTA */

      implicit none

      integer, intent(in)  :: istep

      integer            :: iscf
      logical            :: last, SCFconverged, ts_will_continue
      type(converger_t)  :: conv_harris, conv_etot, conv_free_energy
      real(dp)           :: E_Harris, E_KS_Good
      real(dp)           :: dDmax, dHmax
      real(dp), allocatable :: Hsave(:,:), fa_old(:,:)

#ifdef TRANSIESTA
! Used in force constant runs
      ts_istep = istep
#endif /* TRANSIESTA */
      
!---------------------------------------------------------------------- BEGIN

!     Initialization tasks for a given geometry, including
!     the determination of DM_in
      call state_init( istep )

!     The current structure of the loop tries to reproduce the
!     historical Siesta usage. It should be made more clear.

!     Start of SCF loop
      iscf  = 1
      SCFconverged = .false.

      if (require_free_energy_convergence) then
        call reset(conv_free_energy)
        call set_tolerance(conv_free_energy,freeEnergy_tolerance)
      endif
      if (require_harris_convergence) then
        call reset(conv_harris)
        call set_tolerance(conv_harris,harris_tolerance)
      endif
      if (require_energy_convergence) then
        call reset(conv_etot)
        call set_tolerance(conv_etot,energy_tolerance)
      endif

!     The SCF loop is executed at least once, except if nscf <= 0.
!     Trap this for now. It is not clear what the use-case is...

      if (nscf <= 0)  call die(" *** nscf <= 0.... ")

      call setup_H0()

      ! Beginning of SCF loop

      call timer( 'IterSCF', 1 )
      if (cml_p)
     &     call cmlStartStep( xf=mainXML, type='SCF', index=iscf )

      ! Note that, if mixing H, H_in could also be read from file ...
      H(:,:) = 0.0_dp
      call setup_hamiltonian( iscf ) ! Gets E_KS(DM_in) (non-variational) and H_in
      call compute_max_diff(Hold,H,dHmax)

      DO    ! conditions for exit inside

        call compute_dm( iscf ) ! Gets DM_out (Dscf) and keeps DM_in (Dold)
        call compute_max_diff(Dold,Dscf,dDmax)
        call compute_EKS_Harris( E_Harris, E_KS_Good )

              ! Note: For DM and H convergence checks. At this point: 
              ! If mixing the DM:
              !        Dscf=DM_out, Dold=DM_in, H=H_in, Hold=H_in(prev step)
              !        dDmax=maxdiff(DM_out,DM_in)
              !        dHmax=maxdiff(H_in - H_in(prev step))
              ! If mixing the Hamiltonian:
              !        Dscf=DM_out, Dold=DM_in, H=H_in (mixed in general), Hold=H_in 
              !        dDmax=maxdiff(DM_out,DM_in)
              !        dHmax=maxdiff(H_in - H(DM_out)) ** in prev step **

        call check_scf_convergence( iscf, SCFconverged,
     $                              E_Harris, E_KS_Good,
     $                              dDmax, dHmax,
     $                              conv_harris, conv_etot,
     $                              conv_free_energy)
        if (SCFconverged) then
           if (ionode) then
              write(6,"(a,i0,a)") "SCF cycle converged in ", iscf,
     $                             " iterations."
           endif
        endif

        if (monitor_forces_in_scf) then
           allocate(Hsave(size(H,dim=1),size(H,dim=2)))
           allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
           Hsave = H           ! in case we are mixing H
           fa_old = fa
           call final_H_f_stress(iscf)
           H = Hsave           ! restore
           if (ionode) then
              print *, "Max diff in force (eV/Ang): ",
     $                 maxval(abs(fa-fa_old))*Ang/eV
           endif
           deallocate(Hsave,fa_old)
           call siesta_write_forces()
           call siesta_write_stress_pressure()
        endif

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. (iscf==1)) then
           call timer( 'IterSCF', 3 )
        endif
        if (cml_p) call cmlEndStep(mainXML)

        ! Now check whether we should exit the loop early
        ! In TranSiesta runs, SCF convergence is the signal to switch
        ! to GF mode. Until we completely integrate TranSiesta, we
        ! need to check this.

        if (iscf == nscf)  EXIT
        if (SCFconverged) then
           call check_transiesta_switch(ts_will_continue)
           if (ts_will_continue) then
              ! keep in cycle
           else
              EXIT
              ! Note: we exit with:
              ! If mixing the DM:
              !        Dscf=DM_out, Dold=DM_in, H=H_in, Hold=H_in(prev step)
              ! If mixing the Hamiltonian:
              !        Dscf=DM_out, Dold=DM_in, H=H_in, Hold=H_in 
           endif
        endif

        ! This is formally part of the next SCF step
        ! only done if we are not converged


        call timer( 'IterSCF', 1 )
        if (cml_p) call cmlStartStep(mainXML)

        if (mixH) then
           call setup_hamiltonian( iscf ) ! Gets H_out (H) and keeps H_in (Hold)
           call compute_max_diff(Hold,H,dHmax)
           call mixer (iscf)    ! Mixes H_in and H_out to give H_in_next
        else
           call mixer (iscf)    ! Mixes DM_in and DM_out to give DM_in_next
           call setup_hamiltonian (iscf) ! Gets H_in_next from DM_in_next, and E_KS(DM_in)
           call compute_max_diff(Hold,H,dHmax)
        endif

        ! Save mixed DM (and maybe H later) here for restarts
        call save_density_matrix("mixed")

        iscf = iscf + 1

      ENDDO

      !     End of standard SCF loop.

      if ((SCFMustConverge).and.(.not.SCFconverged)) then
        call die('SCF did not converge in maximum number of steps.')
      endif


      ! Save DM for further analysis offline
      call save_density_matrix("dmout")

      ! Do one more pass to compute forces and stresses (and energies)
      ! with DM_out coming out of the cycle
      call post_scf_work( iscf )

      ! Perform any other analyses of the electronic structure
      ! for this geometry
      call state_analysis( istep )

!------------------------------------------------------------------------ END
      END subroutine siesta_forces

      subroutine check_transiesta_switch(ts_continues_GF)

#ifdef TRANSIESTA
      use m_ts_global_vars,      only: TSinit, TSrun, TSiscf
      use m_ts_options,          only: TSmode
      use m_pulay,               only: init_pulay_arrays
      use parallel,              only: IONode
#endif /* TRANSIESTA */

      logical, intent(out) :: ts_continues_GF

#ifndef TRANSIESTA
      ! We are not compiled as TranSiesta
      ! Nothing else to do in the loop.
      ts_continues_GF = .false.
#else
      ts_continues_GF = .false.

      if (TSmode .and. TSinit) then

         ! We're  done with the initial diagon run
         ! Now we start the TRANSIESTA (Green's Functions) run

         TSinit = .false.       
         TSrun = .true.
         TSiscf = 0

         ts_continues_GF = .true.

           !
           ! From now on, a new mixing cycle starts,
           ! so we have to reset the Pulay record keeping
           !
         call init_pulay_arrays()  
           !
         if(IONode) then
           write(6,'(/a)')
     &'                          ************************       '
           write(6,'(a)')
     &'                          *   TRANSIESTA BEGIN   *       '
           write(6,'(a)')
     &'                          ************************       '
         end if
      end if

#endif /* TRANSIESTA */
      end subroutine check_transiesta_switch

      END module m_siesta_forces
