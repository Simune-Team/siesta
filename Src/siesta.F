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
      implicit none

      integer :: istep, iscf
      logical :: first, last, relaxd

      external:: automatic_cell ! Finds automatic unit cell
      external:: bonds          ! Finds interatomic "bonds" (close atom pairs)
      external:: cgvc           ! CG geometry relaxation including unit cell
      external:: cgvc_zmatrix   ! CG geometry relaxation using Z-matrix coordinates
      external:: fixed          ! Orthogonalizes atomic forces & stress to constraints
      external:: dhscf          ! Finds SCF-dependent part of the Hamiltonian matrix
      external:: diagon         ! Solves the generalized eigenvalue problem
      external:: dnaefs         ! Finds a correction to the neutral-atom energy&forces
      external:: eggbox         ! Subtracts eggbox effect from energy and forces
      external:: extrapol       ! Extrapolates DM for new geometry
      external:: initatom       ! Initializes pseudoatom info
      external:: ioeig          ! I/O of electron eigenenergies
      external:: iofa           ! I/O of atomic forces
      external:: iomd           ! I/O of atomic positions and velocities in MD
      external:: iokp
      external:: iodm           ! I/O of density matrix
      external:: iozm           ! I/O of Z-matrix
      external:: kinefsm        ! Finds kinetic-energy part of the Hamiltonian
      external:: memory         ! Counts allocated memory sizes
      external:: mulliken       ! Finds Mulliken population charges
      external:: naefs          ! Finds neutral-atom part of Hamiltonian, energy&forces
      external:: neighb         ! Finds lists of neighbor atoms
      external:: pulayx         ! Mixes density matrix using Pulay algorithm
      external:: reinit         ! Initializes data read
      external:: shaper         ! Finds geometry type
      external:: spnvec         ! Finds spin vector
      external:: timer          ! Measures and prints CPU time
      external:: xijorb         ! Finds vectors between basis-orbital centers
      external:: prversion


!------------------------------------------------------------------------- BEGIN
      call siesta_init()

      ! Begin of coordinate relaxation iteration
      ! Notice that this loop is not indented
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
