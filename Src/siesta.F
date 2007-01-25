      Program SIESTA
      USE SIESTA_TODO
      implicit none

      integer :: istep, iscf
      logical :: first, last, relaxd

      external :: automatic_cell,
     .  bonds, cgvc, cgvc_zmatrix, fixed,
     .  dhscf, diagon, dnaefs, extrapol, initatom,
     .  iodm, iozm,
     .  kgrid, kgridinit, kinefsm, mulliken, naefs, neighb,
     .  pulayx, 
     .  reinit, shaper, spnvec, 
     .  timer, xijorb, memory,
     .  ioeig, iofa, iokp, iomd, prversion, eggbox


!------------------------------------------------------------------------- BEGIN
      call siesta_init()

      ! Start loop over coordinate changes
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

!       Logica alternativa si no queremos hacer el bucle SCF ni una sola vez
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

        call siesta_move( istep, relaxd )
        istep = istep + 1
      ENDDO
      ! End of coordinate-relaxation loop 

      call siesta_analysis( first, relaxd, iscf )

      call siesta_end()
!--------------------------------------------------------------------------- END
      END program siesta
