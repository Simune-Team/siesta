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
#ifdef MPI      
      use mpi_siesta
#endif      
      use precision, only: dp
      use files, only: slabel
      use siesta_cml
#ifdef SIESTA__FLOOK
      use flook_siesta, only : slua_call, LUA_SCF_LOOP
      use siesta_dicts
#endif
      use m_state_init
      use m_setup_hamiltonian
      use m_setup_H0
      use m_compute_dm
      use m_compute_max_diff
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer, only: mixer
      use m_mixing_scf, only: mixing_scf_converged
      use m_mixing_scf, only: mixers_scf_history_init

      use m_rhog,                only: mix_rhog, compute_charge_diff
      use siesta_options
      use parallel,     only : IOnode, SIESTA_worker
      use m_state_analysis
      use m_steps
      use sys, only : die, bye
      use sparse_matrices, only: H, Hold, Dold, Dscf, Eold, Escf
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance
      use siesta_geom,   only: na_u           ! Number of atoms in unit cell
      use m_energies,    only: Etot           ! Total energy
      use m_forces,      only: fa, cfa        ! Forces and constrained forces
      use m_stress,      only: cstress        ! Constrained stress tensor
      use siesta_master, only: forcesToMaster ! Send forces to master prog
      use siesta_master, only: siesta_server  ! Is siesta a server?
      use m_save_density_matrix, only: save_density_matrix
      use m_iodm_old,    only: write_spmatrix
!
      use units,         only: eV, Ang
      use m_pexsi_solver,        only: prevDmax
      use m_forces,              only: fa
      use units,                 only: eV, Ang
      use write_subs,            only: siesta_write_forces
      use write_subs,            only: siesta_write_stress_pressure
#ifdef NCDF_4
      use dictionary
      use m_ncdf_siesta, only : cdf_init_file, cdf_save_settings
      use m_ncdf_siesta, only : cdf_save_state, cdf_save_basis
#endif
      use m_compute_energies, only: compute_energies
      use m_mpi_utils, only: broadcast
      use fdf
#ifdef SIESTA__PEXSI      
      use m_pexsi, only: pexsi_finalize_scfloop
#endif     

#ifdef TRANSIESTA
      use m_ts_options, only : N_Elec
      use m_ts_method
      use m_ts_global_vars,      only: TSmode, TSinit, TSrun
      use siesta_geom,           only: nsc, na_u, xa, ucell, isc_off
      use atomlist,              only: no_u, lasto, Qtot
      use sparse_matrices,       only: sparse_pattern, block_dist
      use sparse_matrices,       only: Escf, S, maxnh
      use m_spin,                only: nspin
      use m_ts_charge, only : ts_get_charges
      use m_ts_charge,           only: TS_RHOCORR_METHOD
      use m_ts_charge,           only: TS_RHOCORR_FERMI
      use m_ts_charge,           only: TS_RHOCORR_FERMI_TOLERANCE
      use m_transiesta,          only: transiesta
      use kpoint_grid, only : gamma_scf
      use m_energies, only : Ef
#endif /* TRANSIESTA */

      implicit none

      integer, intent(in)  :: istep
      integer :: iscf
      logical :: first, last_step
      logical :: SCFconverged
      real(dp) :: dDmax     ! Max. change in DM elements
      real(dp) :: dHmax     ! Max. change in H elements
      real(dp) :: dEmax     ! Max. change in EDM elements
      real(dp) :: drhog     ! Max. change in rho(G) (experimental)
      type(converger_t)  ::  conv_harris, conv_freeE
!---------------------------------------------------------------------- BEGIN
#ifdef TRANSIESTA
      real(dp) :: Qcur
#endif
#ifdef NCDF_4
      type(dict) :: d_sav
#endif
#ifdef MPI
      integer :: MPIerror
#endif

#ifdef DEBUG
      call write_debug( '    PRE siesta_forces' )
#endif

#ifdef SIESTA__PEXSI      
      ! Broadcast relevant things for program logic
      ! These were set in read_options, called only by "SIESTA_workers".
      call broadcast(nscf,comm=true_MPI_Comm_World)
#endif
      
      if (SIESTA_worker)  then
      !     Initialization tasks for a given geometry
         call state_init( istep )
      endif

#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after state_init")
#endif

      if (fdf_get("Sonly",.false.)) then
         if (SIESTA_worker) then
            call timer( 'all', 2 )
            call timer( 'all', 3 )
         endif
         call bye("S only")
      endif

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

#ifdef SIESTA__FLOOK
      ! Add the iscf constant to the list of variables
      ! that are available only in this part of the
      ! routine.
      call dict_variable_add('SCF.iteration',iscf)
      call dict_variable_add('SCF.dD',dDmax)
      call dict_variable_add('SCF.dH',dHmax)
      call dict_variable_add('SCF.dE',dEmax)
      call dict_variable_add('SCF.drhoG',drhog)
#endif

!     Start of SCF loop
      iscf = 0
      
      ! This call computes the non-scf part of H and initializes the
      ! real-space grid structures.  It might be better to split the two,
      ! putting the grid initialization into state_init and moving the
      ! calculation of H_0 to the body of the loop, done if first=.true.  This
      ! would suit "analysis" runs in which nscf = 0
      if (SIESTA_worker) call setup_H0()

#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after setup_H0")
#endif

#ifdef NCDF_4
      ! Initialize the NC file
      if ( write_cdf ) then

         ! Initialize the file...
         call cdf_init_file(trim(slabel)//'.nc',is_MD = .false.)
#ifdef MPI
         call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

         ! Save the settings
         call cdf_save_settings(trim(slabel)//'.nc')
#ifdef MPI
         call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

         d_sav = ('sp'.kv.1)//('S'.kv.1)
         d_sav = d_sav//('nsc'.kv.1)//('xij'.kv.1)
         d_sav = d_sav//('xa'.kv.1)//('cell'.kv.1)
         d_sav = d_sav//('isc_off'.kv.1)
         call cdf_save_state(trim(slabel)//'.nc',d_sav)
         call delete(d_sav)

         ! Save the basis set
         call cdf_save_basis(trim(slabel)//'.nc')

      end if
#endif

      ! The dHmax variable only has meaning for Hamiltonian
      ! mixing, or when requiring the Hamiltonian to be converged.
      dHmax = -1._dp
      dEmax = -1._dp

      DO 
        ! Conditions of exit:
        !  -- At the top, to catch a non-positive nscf
        !  -- At the bottom, based on convergence or # of iterations

#ifdef SIESTA__FLOOK
        ! Communicate with lua
        call slua_call(LUA, LUA_SCF_LOOP)
#endif

        iscf = iscf + 1
        if ( iscf > nscf ) EXIT

        ! Note implications for TranSiesta when mixing H
        ! Now H will be recomputed instead of simply being
        ! inherited, however, this is required as the 
        ! if we have bias calculations as the electric
        ! field across the junction needs to be present.
        first = (iscf == 1)

        if (SIESTA_worker) then
           if (first) then
              if ( converge_Eharr ) then
                 call reset(conv_harris)
                 call set_tolerance(conv_harris,tolerance_Eharr)
              endif
              if ( converge_FreeE ) then
                 call reset(conv_FreeE)
                 call set_tolerance(conv_FreeE,tolerance_FreeE)
              endif
           endif

          call timer( 'IterSCF', 1 )
          if (cml_p)
     &         call cmlStartStep( xf=mainXML, type='SCF', index=iscf )

       endif ! SIESTA_worker

       if (mixH) then
          if ( first .and. SIESTA_worker) then
             if (fdf_get("Read-H-from-file",.false.)) then
                call get_H_from_file()
             else
                call setup_hamiltonian( iscf )
             endif
          endif
          
          call compute_dm( iscf )
          if (SIESTA_worker) then
             ! Maybe set Dold to zero if reading charge or H...
             call compute_max_diff(Dold,Dscf,dDmax)
             if ( converge_EDM ) 
     &            call compute_max_diff(Eold,Escf,dEmax)
             call setup_hamiltonian( iscf )
             call compute_max_diff(Hold,H,dHmax)
          end if
       else
          if (SIESTA_worker) then
             call setup_hamiltonian( iscf )
             call compute_max_diff(Hold,H,dHmax)
          end if
          call compute_dm( iscf )
          if (SIESTA_worker) then
             call compute_max_diff(Dold,Dscf,dDmax)
             if ( converge_EDM ) 
     &            call compute_max_diff(Eold,Escf,dEmax)
          end if
       endif


       if (SIESTA_worker) then

          call compute_energies(iscf)
          if (mix_charge) then
             call compute_charge_diff(drhog)
          endif
        ! Note: For DM and H convergence checks. At this point: 
        ! If mixing the DM:
        !        Dscf=DM_out, Dold=DM_in(mixed), H=H_in, Hold=H_in(prev step)
        !        dDmax=maxdiff(DM_out,DM_in)
        !        dHmax=maxdiff(H_in - H_in(prev step))
        ! If mixing the Hamiltonian:
        !        Dscf=DM_out, Dold=DM_in, H=H_(DM_out), Hold=H_in(mixed)
        !        dDmax=maxdiff(DM_out,DM_in)
        !        dHmax=maxdiff(H(DM_out),H_in)
        call scfconvergence_test( first, iscf, 
     &                            dDmax, dHmax, dEmax,
     &                            conv_harris, conv_freeE,
     &                            SCFconverged )
        !
        ! ** Check this heuristic
        if ( MixH ) then
           prevDmax = dHmax
        else
           prevDmax = dDmax
        endif

        if (SCFconverged .and. (iscf < min_nscf)) then
           SCFconverged = .false.
           if ( IONode ) then
              write(6,"(a,i0)")
     $         "SCF cycle continued for minimum number of iterations: ",
     $            min_nscf
           end if
        end if

        ! Check whether we should step to the next
        ! mixer
        call mixing_scf_converged( SCFconverged )
        
#ifdef TRANSIESTA
        ! In case the user has requested a Fermi-level correction
        ! Then we start by correcting the fermi-level
        if ( SCFconverged .and. TSrun .and.
     &       TS_RHOCORR_METHOD == TS_RHOCORR_FERMI ) then

           call ts_get_charges(N_Elec, block_dist, sparse_pattern,
     &          nspin, maxnh, Dscf, S, Qtot = Qcur )
           
           if ( abs(Qcur - Qtot) > TS_RHOCORR_FERMI_TOLERANCE ) then

              ! Call transiesta with fermi-correct
              call transiesta(iscf,nspin, 
     &             block_dist, sparse_pattern, Gamma_Scf, ucell, nsc,
     &             isc_off, no_u, na_u, lasto, xa, maxnh, H, S,
     &             Dscf, Escf, Ef, Qtot, .true.) 
           
              ! We will not have not converged as we have just
              ! changed the Fermi-level
              SCFconverged = .false.

           end if
           
        end if
#endif

        if (monitor_forces_in_scf) call compute_forces()

        ! Mix_after_convergence preserves the old behavior of
        ! the program.
        if ((.not. SCFconverged) .or. mix_after_convergence) then
           ! Mix for next step
           if ( mix_charge ) then
              call mix_rhog(iscf)
           else
              call mixer( iscf )
           endif
           ! Save for possible restarts
           if (MixH) then
              call write_spmatrix(H,file="H_MIXED",when=writeH)
              call save_density_matrix(file="DM_OUT",when=writedm)
           else
              call save_density_matrix(file="DM_MIXED",when=writedm)
              call write_spmatrix(H,file="H_DMGEN",when=writeH)
           endif
        endif

        call timer( 'IterSCF', 2 )
        call maybe_print_timings()
        if (cml_p) call cmlEndStep(mainXML)

        ! This should be the last step if we have converged...
        last_step = SCFConverged

#ifdef TRANSIESTA
        ! ... except that we might continue for TranSiesta
        if ( last_step ) then
           call transiesta_switch() ! might reset 'last_step'
                                    ! and iscf
        end if
#endif

        if ( iscf == nscf ) then
           last_step = .true.
        end if

       endif ! SIESTA_worker

#ifdef SIESTA__PEXSI       
       call broadcast(last_step,comm=true_MPI_Comm_World)
       call broadcast(iscf,comm=true_MPI_Comm_World)
#endif
       if ( last_step ) exit

      end do

#ifdef SIESTA__PEXSI      
      if (isolve == SOLVE_PEXSI) then
         call pexsi_finalize_scfloop()
      endif
#endif            

      if (.not. SIESTA_worker) RETURN

      call end_of_cycle_save_operations()
      
      if (.not.SCFconverged) then
         if (ionode.and..not.harrisfun)
     &        write(6,"(a)") ':!: SCF did not converge'//
     $                          ' in maximum number of steps.'
         if (SCFMustConverge) call die('SCF did not converge' //
     $                          ' in maximum number of steps.')
         
      end if

#ifdef TRANSIESTA
      if ( TSmode.and.TSinit.and.(.not. SCFConverged) ) then
         call die('SCF did not converge before proceeding to transiesta
     &calculation')
      end if
#endif

      ! Clean-up here to limit memory usage
      call mixers_scf_history_init( )

      
!     End of standard SCF loop.
!     Do one more pass to compute forces and stresses

!     Note that this call will no longer overwrite H while computing the
!     final energies, forces and stresses...

      if (fdf_get("compute-forces",.true.)) then
         call post_scf_work( istep, iscf , SCFconverged )
#ifdef SIESTA__PEXSI
         if (ionode) call memory_snapshot("after post_scf_work")
#endif
      endif

!     ... so H at this point is the latest generator of the DM, except
!     if mixing H beyond self-consistency or terminating the scf loop
!     without convergence while mixing H

      call state_analysis( istep )
#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after state_analysis")
#endif

      ! If siesta is running as a subroutine, send forces to master program
      if (siesta_server) 
     &  call forcesToMaster( na_u, Etot, cfa, cstress )

#ifdef DEBUG
      call write_debug( '    POS siesta_forces' )
#endif

      CONTAINS

!-----------------------------------------------------
      subroutine get_H_from_file()
      use sparse_matrices, only: maxnh, numh, listh, listhptr
      use atomlist,        only: no_l
      use m_spin,          only: nspin
      use m_iodm_old,      only: read_spmatrix

      logical :: found

      call read_spmatrix(maxnh, no_l, nspin, numh,
     .     listhptr, listh, H, found, userfile="H_IN")
      if (.not. found) call die("Could not find H_IN")
      end subroutine get_H_from_file

!-----------------------------------------------------
      subroutine compute_forces()
      ! Computes forces and stresses with the current DM_out

      use siesta_options, only: recompute_H_after_scf
      use m_final_H_f_stress, only: final_H_f_stress
      use write_subs

      real(dp), allocatable  :: fa_old(:,:), Hsave(:,:)

      allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
      fa_old(:,:) = fa(:,:)
      if ( recompute_H_after_scf ) then
         allocate(Hsave(size(H,dim=1),size(H,dim=2)))
         Hsave(:,:) = H(:,:)
      end if
      call final_H_f_stress( istep , iscf , .false. )
      if ( recompute_H_after_scf ) then
         H(:,:) = Hsave(:,:)
         deallocate(Hsave)
      end if
      if (ionode) then
         print *, "Max diff in force (eV/Ang): ",
     $        maxval(abs(fa-fa_old))*Ang/eV
         call siesta_write_forces(-1)
         call siesta_write_stress_pressure()
      endif
      deallocate(fa_old)

      end subroutine compute_forces

      subroutine maybe_print_timings()
#ifdef TRANSIESTA
      use m_ts_global_vars, only : TSrun
#endif

#ifdef TRANSIESTA
      if ( TSrun ) then
         if ( iscf == 1 ) then
            call timer( 'TS', 3 ) ! with Green function generation
         end if
      else if ( first ) then
         if ( istep == inicoor ) call timer( 'IterSCF', 3 )    
      endif
#else
      if (istep == inicoor .and. first ) call timer( 'IterSCF', 3 )    
#endif

      end subroutine maybe_print_timings
      !
      ! Depending on various conditions, save the DMin
      ! or the DMout, and possibly keep a copy of H

      ! NOTE: Only if the scf cycle converged before exit it
      ! is guaranteed that the DM is "pure out" and that
      ! we can recover the right H if mixing H.
      !
      subroutine end_of_cycle_save_operations()

      if (mix_after_convergence) then
            ! If we have been saving them, there is no point in doing
            ! it one more time
         if (MixH) then
            call save_density_matrix(file="DM_OUT",
     $          when=((.not. writedm) .and. write_dm_at_end_of_cycle))
            call write_spmatrix(H,file="H_MIXED", 
     $          when=((.not. writeH) .and. write_H_at_end_of_cycle))
         else
            call save_density_matrix(file="DM_MIXED",
     $          when=((.not. writedm) .and. write_dm_at_end_of_cycle))
            call write_spmatrix(H,file="H_DMGEN", 
     $          when=((.not. writeH) .and. write_H_at_end_of_cycle))
         endif
      else
         call save_density_matrix(file="DM_OUT",
     $                            when=write_dm_at_end_of_cycle)
         if (MixH) then
             ! Make sure that we keep the H actually used 
             ! to generate the last DM, if needed.
             H = Hold
         endif
         call write_spmatrix(H,file="H_DMGEN", 
     $                       when=write_H_at_end_of_cycle)
      endif

      end subroutine end_of_cycle_save_operations

!-----------------------------------------------------
#ifdef TRANSIESTA
      subroutine transiesta_switch()

      use precision,             only: dp
      use parallel,              only: IONode
      use class_dSpData2D
      use class_Fstack_dData1D

      use siesta_options,        only: fire_mix, broyden_maxit
      use siesta_options,        only: dDtol, dHtol

      use sparse_matrices, only : DM_2D, EDM_2D
      use atomlist, only: lasto
      use siesta_geom, only: nsc, isc_off, na_u, xa, ucell
      use m_energies, only : Ef
      use m_mixing_scf, only: scf_mixs, scf_mix
      use m_mixing, only: mixers_history_init

      use m_ts_global_vars,      only: TSinit, TSrun
      use m_ts_global_vars,      only: ts_print_transiesta
      use m_ts_method
      use m_ts_options,          only: N_Elec, Elecs
      use m_ts_options,          only: DM_bulk
      use m_ts_options,          only: val_swap
      use m_ts_options,          only: ts_scf_mixs
      use m_ts_options,          only: ts_Dtol, ts_Htol
      use m_ts_options,          only: ts_hist_keep
      use m_ts_options,          only: ts_siesta_stop
      use m_ts_electype

      integer :: iEl, na_a
      integer, allocatable :: allowed_a(:)
      real(dp), pointer :: DM(:,:), EDM(:,:)

      ! We are done with the initial diagon run
      ! Now we start the TRANSIESTA (Green functions) run
      if ( .not. TSmode ) return
      if ( .not. TSinit ) return

      ! whether we are in siesta initialization step
      TSinit = .false.
      ! whether transiesta is running
      TSrun  = .true.

      ! If transiesta should stop immediately
      if ( ts_siesta_stop ) then
         
         if ( IONode ) then
          write(*,'(a)') 'ts: Stopping transiesta (user option)!'
         end if
         
         return
         
      end if

      ! Signal to continue...
      last_step = .false.
      iscf   = 0

      ! DANGER (when/if going back to the DIAGON run, we should
      ! re-instantiate the original mixing value)
      call val_swap(dDtol,ts_Dtol)
      call val_swap(dHtol,ts_Htol)
      
      ! From now on, a new mixing cycle starts,
      ! Check in mixer.F for new mixing schemes.
      ! NP new mixing
      if ( associated(ts_scf_mixs, target=scf_mixs) ) then
         do iel = 1 , size(scf_mix%stack)
            call reset(scf_mix%stack(iel), -ts_hist_keep)
            ! Reset iteration count as certain 
            ! mixing schemes require this for consistency
            scf_mix%cur_itt = n_items(scf_mix%stack(iel))
         end do
      else
         call mixers_history_init(scf_mixs)
      end if
      ! Transfer scf_mixing to the transiesta mixing routine
      scf_mix => ts_scf_mixs(1)
      
      call ts_print_transiesta()

      ! In case of transiesta and DM_bulk.
      ! In case we ask for initialization of the DM in bulk
      ! we read in the DM files from the electrodes and 
      ! initialize the bulk to those values
      if ( DM_bulk > 0 ) then
         
         if ( IONode ) then
            write(*,'(/,2a)') 'transiesta: ',
     &           'Initializing bulk DM in electrodes.'
         end if
         
         na_a = 0
         do iEl = 1 , na_u
            if ( .not. a_isDev(iEl) ) na_a = na_a + 1
         end do
         allocate(allowed_a(na_a))
         na_a = 0 
         do iEl = 1 , na_u
            ! We allow the buffer atoms as well (this will even out the
            ! potential around the back of the electrode)
            if ( .not. a_isDev(iEl) ) then
               na_a = na_a + 1
               allowed_a(na_a) = iEl
            end if
         end do
          
         do iEl = 1 , N_Elec
             
            if ( IONode ) then
               write(*,'(/,2a)') 'transiesta: ',
     &              'Reading in electrode TSDE for '//
     &              trim(Elecs(iEl)%Name)
            end if

            ! Copy over the DM in the lead
            ! Notice that the EDM matrix that is copied over
            ! will be equivalent at Ef == 0
            call copy_DM(Elecs(iEl),na_u,xa,lasto,nsc,isc_off, 
     &           ucell, DM_2D, EDM_2D, na_a, allowed_a)

         end do

         ! Clean-up
         deallocate(allowed_a)

         if ( IONode ) then
            write(*,*) ! new-line
         end if

         ! The electrode EDM is aligned at Ef == 0
         ! We need to align the energy matrix
         DM  => val(DM_2D)
         EDM => val(EDM_2D)
         iEl =  size(DM)
         call daxpy(iEl,Ef,DM(1,1),1,EDM(1,1),1)

      end if

      end subroutine transiesta_switch
#endif

      END subroutine siesta_forces

      END module m_siesta_forces
