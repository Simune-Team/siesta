      module m_siesta_forces
      private
      
      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      use siesta_cml
      use m_state_init
      use m_setup_hamiltonian
      use m_compute_dm
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use parallel, only: IOnode
      use m_state_analysis
      use m_steps
      use sys, only : die
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance
#ifdef TRANSIESTA
C TSS Begin
C SIESTA Modules
      use sparse_matrices, only : nh,Dscf, H
      use m_spin, only : nspin
      use sys, only : die
C  TS Modules
      use m_ts_options, only : VIn, VOut, mixH, initH
C TSS End
#endif /* TRANSIESTA */


      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: first, last, SCFconverged
      real(dp) :: dDmax     ! Max. change in density matrix elements
      type(converger_t)  ::  conv_harris, conv_etot
#ifdef TRANSIESTA
      ! TSS Begin
      ! Internal Variables
      integer :: ts_alloc_stat
      ! TSS End
#endif /* TRANSIESTA */
!------------------------------------------------------------------------- BEGIN

        ! Initialization tasks for a given geometry
        call state_init( istep )

        ! The current structure of the loop tries to reproduce the
        ! historical Siesta usage. It should be made more clear.

        ! Start of SCF loop

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

        ! The SCF loop is executed at least once, 
        ! except if nscf <= 0.
        if (nscf <= 0) then
           last = .true.
        endif
#ifdef TRANSIESTA
! TSS Begin
        initH=.true.
! TSS End
#endif /* TRANSIESTA */

        DO WHILE (.not. last)
          call timer( 'IterSCF', 1 )

#ifndef TRANSIESTA
          if (cml_p)
     $         call cmlStartStep(xf=mainXML, type='SCF', index=iscf)

          call setup_hamiltonian( first, last, iscf )

          call compute_dm( iscf )

          call mixer( iscf, dDmax )
#else /* TRANSIESTA */
! TSS Begin: Mixing H or DM
          if(size(VIn) /= nh ) then
            if (allocated(VIn)) then
              deallocate(VIn)
            end if
            if (allocated(VOut)) then
              deallocate(VOut)
            end if

            allocate(VIn(nh,nspin), STAT=ts_alloc_stat)
            if(ts_alloc_stat /= 0)
     . call die('Failled to Allocate VIn in siesta_forces !!')
            allocate(VOut(nh,nspin))
            if(ts_alloc_stat /= 0)
     . call die('Failled to Allocate VOut in siesta_forces !!')

          endif
! TSS End

! Adding the possibility of Mixinh H
          if (mixH) then

           if ( first ) then
            call setup_hamiltonian( first, last, iscf )
            initH=.false.
           end if

           if (cml_p)
     $         call cmlStartStep(xf=mainXML, type='SCF', index=iscf)

           VIn(1:nh,1:nspin)=H(1:nh,1:nspin)

           call compute_dm( iscf )

           call setup_hamiltonian( first, last, iscf )

           VOut(1:nh,1:nspin)=H(1:nh,1:nspin)

           call mixer( iscf, dDmax )

           H(1:nh,1:nspin)=VIn(1:nh,1:nspin)

          else
! TSS End Mixing H

! TSS Begin
! Putting VIn=rho
           VIn(1:nh,1:nspin)=Dscf(1:nh,1:nspin)
           initH=.false.
! TSS End
           if (cml_p)
     $         call cmlStartStep(xf=mainXML, type='SCF', index=iscf)

           call setup_hamiltonian( first, last, iscf )

           call compute_dm( iscf )

! TSS Begin
! Putting VOut=rho
           VOut(1:nh,1:nspin)=Dscf(1:nh,1:nspin)
! TSS End


           call mixer( iscf, dDmax )

! TSS Begin
! Putting rho=VIn after mixing
           Dscf(1:nh,1:nspin)=VIn(1:nh,1:nspin)
! TSS End


! TSS Begin
         end if ! mixH
! TSS End
#endif /* TRANSIESTA */

          call scfconvergence_test( first, iscf, last, dDmax,
     $                              conv_harris, conv_etot )
          ! At this point last is true if SCF has converged. But this
          ! changed below
          SCFconverged = last

          call timer( 'IterSCF', 2 )
          if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )

          iscf = iscf + 1
          first = .false.

          ! The last scheduled iteration is used to compute forces
          ! and stresses, hence the exit if iscf==nscf.
          ! The ">" clause catches the nscf==1 case. For this,

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

        ! End of standard SCF loop.
        ! Do one more pass to compute forces and stresses

        if (last) then
          call post_scf_work( first, last, iscf )
        else
          STOP 'Bad logic in SCF loop'
        endif

        call state_analysis( istep )

!--------------------------------------------------------------------------- END
        END subroutine siesta_forces
      END module m_siesta_forces
