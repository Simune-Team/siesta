      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      use siesta_cml
      use m_state_init
      use m_setup_hamiltonian
      use m_setup_H0
      use m_compute_dm
      use m_compute_max_diff
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer
      use m_rhog,                only: mix_rhog, compute_charge_diff
      USE siesta_options
      use parallel,     only : IOnode, Nodes
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
     
      use atomlist,      only: no_s, no_l, no_u, qtot, indxuo
      use m_spin,        only: nspin
      use m_gamma
      use Kpoint_grid,    only: nkpnt, kpoint, kweight
      use m_eo

      use units,         only: eV, Ang
      use sparse_matrices, only: H, Hold, Dold, Dscf, maxnh, numh,
     &                           listh, listhptr, S, xijo
      use m_save_density_matrix, only: save_density_matrix
      use m_compute_energies, only: compute_energies
!tddft
      use alloc             
      use densematrix,      only: psi,Saux
      use files,            only: slabel
      use m_initwf,         only: initwf  
      use wavefunctions,    only: wavef_ms, iowavef 
      use atomlist,         only: indxua, rmaxo, lasto, iphorb
      use sparse_matrices,  only: Escf, H_kin, H_vkb
      use m_iodm,           only: read_dynamic_dm
      use m_evolve,         only: evolve
      use elec_dyn_options
      use m_iotddft,        only: write_tddft
      use m_overfsm,        only: overfsm
      use siesta_geom,      only: na_s, xa, isa, scell
      use m_energies
      use m_final_H_f_stress,     only: final_H_f_stress
! End tddft - Rafi


#ifdef TRANSIESTA
      use m_ts_global_vars, only : ts_istep
#endif /* TRANSIESTA */
      implicit none

      integer, intent(inout)  :: istep

      integer  :: iscf
      logical  :: first, last_step, SCFconverged
      logical  :: writetdfw
      real(dp) :: dmax      ! Max. change in DM or H elements
      real(dp) :: drhog     ! Max. change in rho(G) (experimental)
      type(converger_t)  ::  conv_harris, conv_freeE
      double precision :: EF7
!tddft
      integer              :: ispn, nio, itd, istmove, mod_tded, mod_md
      integer              :: istpp,nhs
      logical              :: rstart, found
      double precision     :: dt_tded
      character*70, external  :: paste
      integer :: kk
      complex(dp) :: varaux,varaux2
! End -tddft Rafi


#ifdef TRANSIESTA
! Used in force constant runs
      ts_istep = istep
#endif /* TRANSIESTA */
      
!---------------------------------------------------------------------- BEGIN


!     Initialization tasks for a given geometry
      call state_init( istep )

!     The current structure of the loop tries to reproduce the
!     historical Siesta usage. It should be made more clear.
!     Two changes: 
!
!       -- The number of scf iterations performed is exactly
!          equal to the number specified (i.e., the "forces"
!          phase is not counted as a final scf step)
!
!       -- At the change to a TranSiesta GF run the variable "first"
!          is implicitly reset to "true".


!     Start of SCF loop
      iscf  = 0
      call setup_H0()


      
      if (td_elec_dyn) then     !tddft


!      call read_dynoptions ()
      nhs=no_l*no_u*2
      call re_alloc(psi,1,nhs,name='psi',routine='chgbasis')
      call re_alloc(Saux,1,nhs,name='Haux',routine='chgbasis') 
      call chgbasis(no_s, nspin, nspin, no_l, maxnh, maxnh,
     &                no_u, numh, listhptr, listh,
     &                numh, listhptr, listh, S,
     &                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &                Dscf, no_u,istep,psi,Saux)
      !call de_alloc(psi,'psi','chgbasis')       
      else 

      DO 
        ! Conditions of exit:
        !  -- At the top, to catch a non-positive nscf
        !  -- At the bottom, based on convergence or # of iterations

        iscf = iscf + 1
        if (iscf > nscf) EXIT

        ! Note implications for TranSiesta when mixing H
        ! Now H will be recomputed instead of simply being
        ! inherited
        first = (iscf == 1)
#ifdef TRANSIESTA
        !!     if (TSrun) first = .false.   ! for old (wrong) behavior 
#endif

        if (first) then
           if (require_harris_convergence) then
              call reset(conv_harris)
              call set_tolerance(conv_harris,harris_tolerance)
           endif
           if (require_energy_convergence) then
              call reset(conv_freeE)
              call set_tolerance(conv_freeE,energy_tolerance)
           endif
        endif

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
        call compute_energies(iscf)
        if (mix_charge) then
           call compute_charge_diff(drhog)
        endif

        call scfconvergence_test( first, iscf, dmax,
     &                            conv_harris, conv_freeE,
     &                            SCFconverged )

        if (SCFconverged .and. (iscf < min_nscf)) then
           SCFconverged = .false.
           if (ionode) then
              write(6,"(a,i0)")
     $         "SCF cycle continued for minimum number of iterations: ",
     $            min_nscf
           endif
        endif

        if (monitor_forces_in_scf)  call compute_forces()

        if ((.not. SCFconverged) .OR. mix_after_convergence) then
           ! Mix for next step
           if (mix_charge) then
              call mix_rhog(iscf)
           else
              call mixer( iscf )
           endif
           ! Save for possible restarts
           ! Could save H in the future
           call save_density_matrix("mixed")
        endif

        call timer( 'IterSCF', 2 )

        call maybe_print_timings()

        ! This should be the last step if we have converged...
        last_step = SCFConverged
#ifdef TRANSIESTA
        ! ... except that we might continue for TranSiesta
        call check_transiesta_continuation()   ! might reset 'last_step'
                                               ! and iscf
#endif

        if (iscf == nscf) then
           last_step = .true.
        endif

        if (cml_p) call cmlEndStep(mainXML)

        if (last_step) EXIT

      ENDDO

      if ((SCFMustConverge).and.(.not.SCFconverged)) then
        call die('SCF did not converge in maximum number of steps.')
      endif

      if (.not. mix_after_convergence) then
         ! Save final DM_out for further analysis offline
         call save_density_matrix("dmout")
         if (MixH) then
             ! We could make sure that we keep the H actually used 
             ! to generate the last DM, if needed.
             !!!    H = Hold
         endif
      endif
             !To write the initial wavefunctions to be used in a consequent TDDFT
             ! run.
      if(writetdwf) then
              istpp = 0
      call initwf(no_s, nspin, nspin, no_l, maxnh, no_u,
     &            numh, listhptr, listh, H, S, qtot,
     &            gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     &            eo, qo, no_u, EF7, istpp,totime)
      end if


!     End of standard SCF loop.
!     Do one more pass to compute forces and stresses

      call post_scf_work( iscf )

      call state_analysis( istep )

      ! If siesta is running as a subroutine, send forces to master program
      if (siesta_server) 
     &  call forcesToMaster( na_u, Etot, cfa, cstress )
  
      end if !tddft

      if (td_elec_dyn) then
      ! Read the initial wavefunctions. 
      
      IF (istep==1) THEN
       allocate(wavef_ms(nkpnt,nspin))
       call iowavef('read',wavef_ms,no_u,nkpnt,nspin, istpp,rstart_time)

        istep=istpp+1
      ! istep=istpp
      END IF 

! To save the TD-wavefunctions and Density matrix to restart
        mod_md   = mod(istep, tdednwrite)

        if(mod_md .eq. 0 .or. istep .eq. fincoor) then
           call iowavef('write',wavef_ms,no_u,nkpnt,nspin, istep,totime)
           call save_density_matrix("dmout")
        end if


      DO itded = 1,ntded  !TDED loop

      dt_tded= dt/ntded
         if(IONode) then
           write(6,'(/a)')
     & '                     ************************       '
           write(6,*)'                 TDED Step    =', itded
           write(6,'(/a)')
     & '                     ************************       '
        end if
 
! To save the TD-wavefunctions and Density matrix to restart
        mod_tded = mod(itded, tdednwrite)

        if(mod_tded .eq. 0) then
          call iowavef('write',wavef_ms,no_u,nkpnt,nspin, istep, totime)
          call save_density_matrix("dmout")
        end if

      call setup_hamiltonian( itded )
      
!      if (istpp == 0) then
      call compute_energies (itded)
      call write_tddft (totime, istep, itded, ntded, rstart_time,      
     &                   etot, eo, no_u,nspin,nkpnt) 
!      end if
           call evolve(no_s,nspin,nspin,no_l,maxnh,maxnh,no_u,
     &            numh,listhptr,listh,numh,listhptr,listh,
     &            H,S,eo,
     &            gamma, xijo, indxuo,nkpnt,kpoint,kweight,
     &            Dscf, Escf, no_u, dt_tded,istep,itded)
!if (IONode)print*,eo
! The total simulation time mainly for plotting
       totime = (istep*dt - dt) + (itded*dt_tded-dt_tded) 

!      call compute_energies (itded)
      call final_H_f_stress (itded)
      call post_scf_work( itded )
      call state_analysis( istep )
      
!      call write_tddft (totime, istep, itded, ntded, rstart_time,      
!     &                   etot, eo, no_u,nspin,nkpnt) 
      
!    if (etot_time) then
!           call ioetot(totime, Etot, istep, itded, ntded,rstart_time)
!           endif
!           if (dip_time) then
!           call write_td_dipol(totime, istep, itded, ntded,rstart_time) 
!           endif
!           if (eigen_time) then
!           call ioeigenvalues (eo,no_u,nspin,nspin,nkpnt, 
!     &                       totime,rstart_time, istp, itded,ntded)
!           endif

!
        END DO !  TDED loop



      end if !    td_elec_dyn





!------------------------------------------------------------------------ END

      CONTAINS

!-----------------------------------------------------
      subroutine compute_forces()
      ! Computes forces and stresses with the current DM_out

      use m_final_H_f_stress,    only: final_H_f_stress
      use write_subs

      real(dp), allocatable  :: fa_old(:,:), Hsave(:,:)

      allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
      allocate(Hsave(size(H,dim=1),size(H,dim=2)))
      fa_old = fa
      Hsave(:,:) = H(:,:)
      call final_H_f_stress(iscf)
      H(:,:) = Hsave(:,:)
      if (ionode) then
         print *, "Max diff in force (eV/Ang): ",
     $        maxval(abs(fa-fa_old))*Ang/eV
         call siesta_write_forces()
         call siesta_write_stress_pressure()
      endif
      deallocate(fa_old,Hsave)

      end subroutine compute_forces

!-----------------------------------------------------
      subroutine maybe_print_timings()
#ifdef TRANSIESTA
      use m_ts_global_vars, only: TSrun, TSiscf
#endif

#ifdef TRANSIESTA
      if (TSrun) then
         if (TSiscf .eq. 1) then
            call timer( 'TS', 3 ) !TSS (with Green's function generation)
            call timer( 'TS_calc', 3 ) !TSS (only calculating routines)
         end if
      else
         if (istep.eq.inicoor .and. first ) call timer( 'IterSCF', 3 )    
      endif
#else
      if (istep.eq.inicoor .and. first ) call timer( 'IterSCF', 3 )    
#endif
      end subroutine maybe_print_timings

!-----------------------------------------------------

#ifdef TRANSIESTA
      subroutine check_transiesta_continuation()

      ! This routine performs some basic control of the
      ! scf loop for TranSiesta calculations. 

      ! Recall: 
      !    TSmode: .true. if we are in TranSiesta mode
      !    TSinit: .true. if performing an initial Siesta phase
      !    TSrun : .true. if in an actual TS GF run.

      use m_ts_global_vars, only: TSinit, TSrun, TSiscf
      use m_ts_options, only: TSmode
      use m_pulay, only: init_pulay_arrays

      ! Signal the change in character of the scf loop

      if (TSmode .and. last_step .and. TSinit) then

         last_step = .false.  ! we're only done with the initial diagon run
         TSinit = .false.     ! now we start the TRANSIESTA run
         TSrun = .true.
         TSiscf = 0
         iscf = 0
         !
         ! From now on, a new mixing cycle starts,
         ! so we have to reset the Pulay record keeping
         !
         call init_pulay_arrays()

         if(IONode) then
           write(6,'(/a)')
     &'                          ************************       '
           write(6,'(a)')
     &'                          *   TRANSIESTA BEGIN   *       '
           write(6,'(a)')
     &'                          ************************       '
        end if
      end if


      ! Avoid stopping the loop too early
      if (TSrun .and. TSiscf < min_nscf ) last_step = .false.
! Needed for an MD run of TranSiesta
!      if (TSrun .and. last_step) then
!         ! reset pulays again
!         call init_pulay_arrays()
!         ! Set TSinit to .true. for a continuation run in an MD simulation.
!         TSinit = .true.
!         !
!         TSrun  = .false.
!         TSiscf = 0
!      end if

      end subroutine check_transiesta_continuation
#endif /* TRANSIESTA */


      END subroutine siesta_forces
      END module m_siesta_forces
