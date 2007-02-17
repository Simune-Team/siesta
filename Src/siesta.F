      Program SIESTA
      USE SIESTA_TODO
      use m_siesta_init
      use m_state_init
      use m_siesta_analysis
      use m_setup_hamiltonian
      use m_siesta_move
      use m_siesta_end
      use m_compute_dm
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use siesta_geom
      use molecularmechanics, only : twobody
      use parallel, only: IOnode
      use m_state_analysis
      use m_stress
      implicit none

      integer :: istep, iscf
      logical :: first, last, relaxd

!------------------------------------------------------------------------- BEGIN
      call siesta_init()

      ! Begin of coordinate relaxation iteration
      relaxd = .false.
      istep = inicoor
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
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

          if (first) then
            ! Add on force field contribution if required
            call twobody(na_u,xa,isa,ucell,Emm,ifa,fa,istr,stress)
          endif

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
        call siesta_move( istep, relaxd )
        if (.not. relaxd) then
          istep = istep + 1
        endif
      ENDDO
      ! End of coordinate-relaxation loop 

      call siesta_analysis( first, relaxd, iscf )

      call siesta_end()
!--------------------------------------------------------------------------- END
      END program siesta
