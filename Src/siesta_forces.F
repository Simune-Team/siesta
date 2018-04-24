! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)

      use precision, only: dp
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
      use parallel,     only : IOnode
      use m_state_analysis
      use m_steps
      use sys, only : bye
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance
      use siesta_geom,   only: na_u           ! Number of atoms in unit cell
      use m_energies,    only: Etot           ! Total energy
      use m_forces,      only: fa, cfa        ! Forces and constrained forces
      use m_stress,      only: cstress        ! Constrained stress tensor
      use siesta_master, only: forcesToMaster ! Send forces to master prog
      use siesta_master, only: siesta_server  ! Is siesta a server?

      use units,         only: eV, Ang
      use sparse_matrices, only: H, Hold, Dold, Dscf
      use m_save_density_matrix, only: save_density_matrix
      use m_compute_energies, only: compute_energies
      use m_check_walltime
      use m_mpi_utils, only: barrier
      
#ifdef TRANSIESTA
      use m_ts_global_vars, only : ts_istep
#endif /* TRANSIESTA */
      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: first, last_step, SCFconverged
      real(dp) :: dmax      ! Max. change in DM or H elements
      real(dp) :: drhog     ! Max. change in rho(G) (experimental)
      type(converger_t)  ::  conv_harris, conv_freeE

      logical               :: time_is_up
      character(len=40)     :: tmp_str
      
      external :: die, message

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

      ! This call computes the non-scf part of H and initializes
      ! the real-space grid structures.
      ! It might be better to split the two, putting the grid initialization
      ! into state_init and moving the calculation
      ! of H_0 to the body of the loop, done if first=.true.
      ! This would suit "analysis" runs in which nscf = 0
      !
      call setup_H0()

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
           SCFConverged = .false.
           if (require_harris_convergence) then
              call reset(conv_harris)
              call set_tolerance(conv_harris,harris_tolerance)
           endif
           if (require_energy_convergence) then
              call reset(conv_freeE)
              call set_tolerance(conv_freeE,energy_tolerance)
           endif
        endif

        ! Check whether we are short of time to continue
        call check_walltime(time_is_up)
        if (time_is_up) then
           ! Save DM/H if we were not saving it...
           !     Do any other bookeeping not done by "die"
           call timer('all',2)
           call timer('all',3)
           if (.not. SCFConverged) then
              call message('WARNING',
     $             'SCF_NOT_CONV: SCF did not converge'//
     $             ' before wall time exhaustion')
              write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, dmax
              call message(' (info)',"Geom step, scf iteration, dmax:"
     $                            //trim(tmp_str))
           endif
           call barrier() ! A non-root node might get first to the 'die' call
           call die("OUT_OF_TIME: Time is up.")
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
     $         ":!: SCF cycle continued for minimum number"//
     $             " of iterations: ", min_nscf
           endif
        endif

        if (monitor_forces_in_scf)  call compute_forces()

        ! Mix_after_convergence preserves the old behavior of
        ! the program.
        if ((.not. SCFconverged) .OR. mix_after_convergence) then
           ! Mix for next step
           if (mix_charge) then
              call mix_rhog(iscf)
           else if (MixH) then
              ! Save the DM (which is actually "out") for possible restarts
              ! In future we might save H
              call save_density_matrix("out_but_not_final")
              call mixer( iscf )
           else 
              call mixer( iscf )
              ! Save for possible restarts
              call save_density_matrix("mixed")
           endif
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

      if (SCFconverged) then

         ! Only if the scf cycle converged before exit it
         ! is guaranteed that the DM is "pure out" and that
         ! we can recover the right H if mixing H.

	 ! There are still more cases that could be handled                                                                
	 ! more transparently with logical state variables
         ! such as h_is_hdmgen, dm_is_dmout, etc...
         ! The PEXSI branch contains more logic.

         if (.not. mix_after_convergence) then
            ! Save final DM_out for further analysis offline
            call save_density_matrix("dmout")
            if (MixH) then
               ! Keep the H actually used to generate the last DM
               H = Hold
            endif
         endif

      else

         if (SCFMustConverge) then
            call message('FATAL','SCF_NOT_CONV: SCF did not converge' //
     $           ' in maximum number of steps (required).')
            write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, dmax
            call message(' (info)',"Geom step, scf iteration, dmax:"
     $           //trim(tmp_str))
            call timer( 'all', 2 ) ! New call to close the tree
            call timer( 'all', 3 )
            call barrier()
            call die('ABNORMAL_TERMINATION')
         else
            call message('WARNING',
     $           'SCF_NOT_CONV: SCF did not converge'//
     $           ' in maximum number of steps.')
            write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, dmax
            call message(' (info)',"Geom step, scf iteration, dmax:"
     $                            //trim(tmp_str))
         endif
      endif


!     End of standard SCF loop.
!     Do one more pass to compute forces and stresses

!     Note that this call will no longer overwrite H while computing the
!     final energies, forces and stresses...
!
      call post_scf_work( iscf )

!     ... so H at this point is the latest generator of the DM, except
!     if mixing H beyond self-consistency or terminating the scf loop
!     without convergence while mixing H
      
      call state_analysis( istep )

      ! If siesta is running as a subroutine, send forces to master program
      if (siesta_server) 
     &  call forcesToMaster( na_u, Etot, cfa, cstress )


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
         call siesta_write_forces(-1)
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
