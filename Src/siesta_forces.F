      module m_siesta_forces
      private
      
      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      use m_state_init
      use m_setup_hamiltonian
      use m_setup_H0
      use m_compute_dm
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use parallel, only: IOnode
      use m_state_analysis
      use m_steps

      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: first, last
      real(dp) :: dDmax     ! Max. change in density matrix elements

!------------------------------------------------------------------------- BEGIN

        ! Initialization tasks for a given geometry
        call state_init( istep )

        ! The current structure of the loop tries to reproduce the
        ! historical Siesta usage. It should be made more clear.

        ! Start of SCF loop

        iscf  = 1
        first = .true.
        last  = .false.

        ! The SCF loop is executed at least once, 
        ! except if nscf <= 0.
        if (nscf <= 0) then
           last = .true.
        endif
        
        call setup_H0()

        DO WHILE (.not. last)
          call timer( 'IterSCF', 1 )

          ! Will use H0 computed above
          call setup_hamiltonian( iscf )

          call timer( 'compute_dm', 1 )
          call compute_dm( iscf )
          call timer( 'compute_dm', 2 )

          call mixer( iscf, dDmax )

          call timer( 'SCFconv', 1 )
          call scfconvergence_test( first, iscf, last, dDmax )
          call timer( 'SCFconv', 2 )

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

        ENDDO

        ! End of standard SCF loop.
        ! Do one more pass to compute forces and stresses

        if (last) then
          call post_scf_work( iscf )
        else
          STOP 'Bad logic in SCF loop'
        endif

        call state_analysis( istep )

!--------------------------------------------------------------------------- END
        END subroutine siesta_forces
      END module m_siesta_forces
