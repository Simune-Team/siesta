      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
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
      use m_mixer
      use m_rhog,                only: mix_rhog, compute_charge_diff
      use siesta_options
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

      use units,         only: eV, Ang
      use sparse_matrices, only: H, Hold, Dold, Dscf
      use m_save_density_matrix, only: save_density_matrix
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
#ifdef TRANSIESTA
      use m_ts_options, only : N_Elec
      use m_ts_method
      use m_ts_global_vars,      only: TSmode, TSinit, TSrun
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
      use m_pulay,           only: init_pulay_arrays
#endif /* TRANSIESTA */
      implicit none

      integer, intent(in)  :: istep

      integer            :: iscf
      logical            :: first, last_step
      logical            :: SCFconverged
      real(dp) :: dDmax     ! Max. change in DM elements
      real(dp) :: dHmax     ! Max. change in H elements
      real(dp) :: drhog     ! Max. change in rho(G) (experimental)
      type(converger_t)  ::  conv_harris, conv_freeE
#ifdef TRANSIESTA
      real(dp) :: Qcur
#endif
#ifdef NCDF_4
      type(dict) :: d_sav
#endif

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

#ifdef SIESTA__FLOOK
      ! Add the iscf constant to the list of variables
      ! that are available only in this part of the
      ! routine.
      call dict_variable_add('SCF.Iteration',iscf)
      call dict_variable_add('SCF.dDmax',dDmax)
      call dict_variable_add('SCF.dHmax',dHmax)
      call dict_variable_add('SCF.drhoG',drhog)
#endif

!     Start of SCF loop
      iscf = 0
      call setup_H0()

#ifdef NCDF_4
      ! Initialize the SIESTA.nc file
      if ( write_cdf ) then

         ! Initialize the file...
         call cdf_init_file('SIESTA.nc',is_MD = .false.)

         ! Save the settings
         call cdf_save_settings('SIESTA.nc')

         d_sav = ('sp'.kv.1)//('S'.kv.1)
         d_sav = d_sav//('nsc'.kv.1)//('xij'.kv.1)
         d_sav = d_sav//('xa'.kv.1)//('cell'.kv.1)
         d_sav = d_sav//('isc_off'.kv.1)
         call cdf_save_state('SIESTA.nc',d_sav)
         call delete(d_sav)

         ! Save the basis set
         call cdf_save_basis('SIESTA.nc')

      end if
#endif

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

        if ( first ) then
           if ( require_harris_convergence ) then
              call reset(conv_harris)
              call set_tolerance(conv_harris,harris_tolerance)
           endif
           if ( require_energy_convergence ) then
              call reset(conv_freeE)
              call set_tolerance(conv_freeE,energy_tolerance)
           endif
        endif

        call timer( 'IterSCF', 1 )
        if (cml_p)
     &         call cmlStartStep( xf=mainXML, type='SCF', index=iscf )
       
        if ( mixH ) then
           if ( first ) call setup_hamiltonian( iscf )
           call compute_dm( iscf )
           ! Maybe set Dold to zero if reading charge or H...
           call compute_max_diff(Dold,Dscf,dDmax) ! Assuming we have Dold
           call setup_hamiltonian( iscf )
           call compute_max_diff(Hold,H,dHmax)
        else 
           call setup_hamiltonian( iscf )
           ! Maybe make sure that Hold = 0 at the beginning?
           if (.not. first) call compute_max_diff(Hold,H,dHmax) 
           call compute_dm( iscf )
           call compute_max_diff(Dold,Dscf,dDmax)
        end if 
        call compute_energies(iscf)
        if ( mix_charge ) then
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

        call scfconvergence_test( first, iscf, dDmax, dHmax, 
     &                            conv_harris, conv_freeE,
     &                            SCFconverged )

        if (SCFconverged .and. (iscf < min_nscf)) then
           SCFconverged = .false.
           if ( IONode ) then
              write(6,"(a,i0)")
     $         "SCF cycle continued for minimum number of iterations: ",
     $            min_nscf
           end if
        end if

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
              call transiesta(iscf,nspin, 
     &             block_dist, sparse_pattern, Gamma_Scf, ucell, nsc,
     &             isc_off, no_u, na_u, lasto, xa, maxnh, H, S,
     &             Dscf, Escf, Ef, Temp, Qtot, .true.) 
           
              ! We will not have not converged as we have just
              ! changed the Fermi-level
              SCFconverged = .false.

           end if
           
        end if
#endif

        if (monitor_forces_in_scf) call compute_forces()

        if ((.not. SCFconverged) .or. mix_after_convergence) then
           ! Mix for next step
           if ( mix_charge ) then
              call mix_rhog(iscf)
           else
              call mixer( iscf )
           end if
           ! Save for possible restarts
           ! Could save H in the future
           call save_density_matrix("mixed")
        end if

        call timer( 'IterSCF', 2 )

        call maybe_print_timings()

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

        if (cml_p) call cmlEndStep(mainXML)

        if ( last_step ) EXIT

      END DO

      if (SCFMustConverge.and.(.not.SCFconverged)) then
        call die('SCF did not converge in maximum number of steps.')
      endif

#ifdef TRANSIESTA
      if ( TSmode.and.TSinit.and.(.not. SCFConverged) ) then
         call die('SCF did not converge before proceeding to transiesta
     &calculation')
      end if
#endif

      if (.not. mix_after_convergence) then
         ! Save final DM_out for further analysis offline
         call save_density_matrix("dmout")
         if ( MixH ) then
             ! We could make sure that we keep the H actually used 
             ! to generate the last DM, if needed.
             !!!    H = Hold
         end if
      end if

!     End of standard SCF loop.
!     Do one more pass to compute forces and stresses

      call post_scf_work( istep , iscf )

      call state_analysis( istep )

      ! If siesta is running as a subroutine, send forces to master program
      if (siesta_server) 
     &  call forcesToMaster( na_u, Etot, cfa, cstress )

      CONTAINS

!-----------------------------------------------------
      subroutine compute_forces()
      ! Computes forces and stresses with the current DM_out

      use m_final_H_f_stress,    only: final_H_f_stress
      use write_subs

      real(dp), allocatable  :: fa_old(:,:), Hsave(:,:)

      allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
      allocate(Hsave(size(H,dim=1),size(H,dim=2)))
!$OMP parallel workshare default(shared)
      fa_old(:,:) = fa(:,:)
      Hsave(:,:) = H(:,:)
!$OMP end parallel workshare
      call final_H_f_stress( istep , iscf , save_HS = .false. )
!$OMP parallel workshare default(shared)
      H(:,:) = Hsave(:,:)
!$OMP end parallel workshare
      if (ionode) then
         print *, "Max diff in force (eV/Ang): ",
     $        maxval(abs(fa-fa_old))*Ang/eV
         call siesta_write_forces()
         call siesta_write_stress_pressure()
      endif
      deallocate(fa_old,Hsave)

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

#ifdef TRANSIESTA
      subroutine transiesta_switch()

      use precision,             only: dp
      use parallel,              only: IONode
      use class_dSpData2D
      use siesta_options,        only: fire_mix, broyden_maxit
      use siesta_options,        only: wmix, dDtol, dHtol

      use sparse_matrices, only : DM_2D, EDM_2D
      use atomlist, only: no_u, lasto
      use siesta_geom, only: nsc, isc_off, na_u, xa, ucell
      use m_pulay,               only: init_pulay_arrays
      use m_energies, only : Ef

      use m_ts_global_vars,      only: TSinit, TSrun
      use m_ts_method
      use m_ts_options,          only: N_Elec, Elecs
      use m_ts_options,          only: DM_bulk
      use m_ts_options,          only: val_swap
      use m_ts_options,          only: ts_wmix, ts_Dtol, ts_Htol
      use m_ts_electype

      integer :: iEl, na_a
      integer, allocatable :: allowed_a(:)
      real(dp), pointer :: DM(:,:), EDM(:,:)

      ! We are done with the initial diagon run
      ! Now we start the TRANSIESTA (Green functions) run
      if ( .not. TSmode ) return
      if ( .not. TSinit ) return

      last_step = .false.
      TSinit = .false.    ! whether we are in siesta initialization step
      TSrun  = .true.     ! whether transiesta is running
      iscf   = 0
      
      ! DANGER (when/if going back to the DIAGON run, we should
      ! re-instantiate the original mixing value)
      call val_swap(wmix,ts_wmix)
      call val_swap(dDtol,ts_Dtol)
      call val_swap(dHtol,ts_Htol)
      
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
     &'                          ************************'
         write(6,'(a)')
     &'                          *   TRANSIESTA BEGIN   *'
         write(6,'(a)')
     &'                          ************************'
      end if

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
