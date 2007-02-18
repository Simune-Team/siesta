      module m_siesta_forces
      private
      
      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      USE SIESTA_TODO
      use m_state_init
      use m_setup_hamiltonian
      use m_compute_dm
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use parallel, only: IOnode
      use m_state_analysis

      implicit none

      integer, intent(in)  :: istep

      integer :: iscf
      logical :: first, last

!------------------------------------------------------------------------- BEGIN

        call state_init( istep )

        ! Start of SCF iteration 
        iscf  = 1
        first = .true.
        last  = .false.
        if (wmix .le. 0.0_dp) then
          if (IOnode) then
            write(6,'(/,a,f15.8)')
     .     'siesta: WARNING: Mixing weight for SCF loop =', wmix
          endif
          last = .true.
        endif

        !Logica alternativa si no queremos hacer el bucle SCF ni una sola vez
        if (nscf .eq. 1) then
          last = .true.
        endif

        DO WHILE (.not. last)
          call timer( 'IterSCF', 1 )

          call setup_hamiltonian( first, last, iscf )

          call compute_dm( iscf )

          call mixer( iscf )

          call scfconvergence_test( first, iscf, last )

          call timer( 'IterSCF', 2 )
          if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )

          iscf = iscf + 1
          first = .false.
          if (iscf .eq. nscf) then
            last = .true.
          endif

        ENDDO
        ! End of SCF iteration

        if (last) then
          call post_scf_work( first, last, iscf )
        else
          STOP 'MALLLL'
        endif

        call state_analysis( istep )
!--------------------------------------------------------------------------- END
        END subroutine siesta_forces
      END module m_siesta_forces
