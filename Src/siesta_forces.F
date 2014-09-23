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
      use siesta_options,        only: SCFMustConverge, nscf, min_nscf
      use siesta_options,        only: mixH
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
      use m_ts_options, only : N_Elec
      use m_ts_method
      use m_ts_global_vars,      only: TSmode, ts_istep, TSrun, TSiscf
      use siesta_geom,           only: nsc, na_u, xa, ucell, isc_off
      use atomlist,              only: no_u, lasto, Qtot
      use sparse_matrices,       only: sparse_pattern, block_dist
      use sparse_matrices,       only: Escf, S, maxnh, xijo
      use m_spin,                only: nspin
      use m_ts_charge, only : ts_get_charges
      use m_ts_charge,           only: TS_RHOCORR_METHOD
      use m_ts_charge,           only: TS_RHOCORR_FERMI
      use m_ts_charge,           only: TS_RHOCORR_FERMI_TOLERANCE
      use m_transiesta,          only: transiesta
      use kpoint_grid, only : gamma_scf
      use m_energies, only : Ef
      use siesta_options, only : temp
      use m_pulay,           only: init_pulay_arrays
#endif /* TRANSIESTA */

      implicit none

      integer, intent(in)  :: istep

      integer            :: iscf
      logical            :: SCFconverged, ts_will_continue
      type(converger_t)  :: conv_harris, conv_etot, conv_free_energy
      real(dp)           :: E_Harris, E_KS_Good
      real(dp)           :: dDmax, dHmax
      real(dp), allocatable :: Hsave(:,:), fa_old(:,:)
#ifdef TRANSIESTA
      real(dp) :: Qcur
#endif

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
#ifdef TRANSIESTA
      TSiscf = 1
#endif
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
        if (SCFconverged .and. (iscf < min_nscf)) then
           SCFconverged = .false.
           if (ionode) then
              write(6,"(a,i0)")
     $         "SCF cycle continued for minimum number of iterations: ",
     $            min_nscf
           endif
        endif

#ifdef TRANSIESTA
        ! In case the user has requested a Fermi-level correction
        ! Then we start by correcting the fermi-level
        if ( TSrun .and. SCFconverged .and.
     &       TS_RHOCORR_METHOD == TS_RHOCORR_FERMI ) then

           call ts_get_charges(N_Elec, block_dist, sparse_pattern,
     &          nspin, maxnh, Dscf, S, Qtot = Qcur )
           
           if ( abs(Qcur - Qtot) > TS_RHOCORR_FERMI_TOLERANCE ) then

              ! First we reset the pulay arrays
              call init_pulay_arrays()
              
              ! Call transiesta with fermi-correct
              call transiesta(TSiscf,nspin, 
     &             block_dist, sparse_pattern, Gamma_Scf, ucell, nsc,
     &             isc_off, no_u, na_u, lasto, xa, maxnh, H, S,
     &             Dscf, Escf, Ef, Temp, Qtot, .true.) 
           
              ! We will not have not converged as we have just
              ! changed the Fermi-level
              SCFconverged = .false.

           end if
           
        end if
#endif

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
#ifdef TRANSIESTA
        if (istep.eq.inicoor .and. (iscf==1.and..not.TSrun)) then
           call timer( 'IterSCF', 3 )
        endif
        if (istep.eq.inicoor .and. (TSiscf==1.and.TSrun)) then
           call timer( 'TS', 3 )
        endif
#else
        if (istep.eq.inicoor .and. (iscf==1)) then
           call timer( 'IterSCF', 3 )
        endif
#endif        
        if (cml_p) call cmlEndStep(mainXML)

        ! Now check whether we should exit the loop early
        ! In TranSiesta runs, SCF convergence is the signal to switch
        ! to GF mode. Until we completely integrate TranSiesta, we
        ! need to check this.

        if (iscf == nscf) EXIT
        ts_will_continue = .false.
        if (SCFconverged) then
           call check_transiesta_switch(ts_will_continue)
           if ( .not. ts_will_continue) then
              if ( mix_after_convergence ) then
                 ! Should we do setup_hamiltonian if mixH?
                 call mixer (iscf)
              end if
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

#ifdef TRANSIESTA
        if ( ts_will_continue ) then
           ! As transiesta is a "new" solution method we should in principle
           ! have two SCF loops. So we will never mix before entering 
           ! transiesta.
           
           call init_transiesta_switch()
           
           ! Initialize the SCF counters
           TSiscf = 0
           iscf = 0

           ! We can now setup the correct hamiltonian for the 
           ! start of the transiesta cycle
           ! I.e. here we also add the bias in the Vscf array
           call setup_hamiltonian( iscf )
           ! Calculate the dHmax for transfering to
           ! the Hamiltonian created with bias
           call compute_max_diff(Hold,H,dHmax)

        end if
#endif

        ! Save mixed DM (and maybe H later) here for restarts
        call save_density_matrix("mixed")
        
        iscf = iscf + 1
#ifdef TRANSIESTA
        if ( TSrun ) TSiscf = TSiscf + 1
#endif

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
      use m_ts_global_vars, only: TSinit, TSmode
#endif

      logical, intent(out) :: ts_continues_GF

#ifndef TRANSIESTA
      ! We are not compiled as TranSiesta
      ! Nothing else to do in the loop.
      ts_continues_GF = .false.
#else
      ts_continues_GF = (TSmode .and. TSinit)
#endif /* TRANSIESTA */

      end subroutine check_transiesta_switch

#ifdef TRANSIESTA
      subroutine init_transiesta_switch()

      use precision,             only: dp
      use siesta_options,        only: fire_mix, wmix, broyden_maxit
      use m_ts_global_vars,      only: TSinit, TSrun
      use m_ts_options,          only: ts_wmix
      use m_pulay,               only: init_pulay_arrays
      use parallel,              only: IONode

      real(dp) :: tmp

      ! We're  done with the initial diagon run
      ! Now we start the TRANSIESTA (Green's Functions) run

      TSinit = .false. 
      TSrun  = .true.
      
      ! DANGER (when going back to the DIAGON run, we should re-instantiate
      ! the original mixing value)
      tmp     = wmix
      wmix    = ts_wmix
      ts_wmix = tmp
      
      ! From now on, a new mixing cycle starts,
      ! Check in mixer.F for new mixing schemes.
      if ( fire_mix ) then
         ! do nothing
      else if ( broyden_maxit == 0 ) then
         ! so we have to reset the Pulay record keeping
         call init_pulay_arrays() 
      else
         ! do nothing
      end if
      ! the Broyden mixing and the Fire-mixing does
      ! not need to initialize, that is done implicitly

      if ( IONode ) then
         write(6,'(/a)')
     &'                          ************************       '
         write(6,'(a)')
     &'                          *   TRANSIESTA BEGIN   *       '
         write(6,'(a)')
     &'                          ************************       '
      end if

      end subroutine init_transiesta_switch
#endif

      END module m_siesta_forces
