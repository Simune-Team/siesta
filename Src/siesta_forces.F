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
      use parallel,     only : IOnode
      use m_state_analysis
      use m_steps
#ifndef RGT
      use sparse_matrices
      use m_spin,   only: nspin
      use atomlist, only: no_l, no_s
      use m_ordern
#endif

      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: first, last
      real(dp) :: dDmax     ! Max. change in density matrix elements
#ifndef RGT
      logical  :: found
#endif
!------------------------------------------------------------------------- BEGIN
#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 10 )
#endif
C     Initialization tasks for a given geometry
#ifdef DEBUG
      call write_debug( '  PRE state_init' )
#endif
      call state_init( istep )

C     The current structure of the loop tries to reproduce the
C     historical Siesta usage. It should be made more clear.

C     Start of SCF loop
      iscf  = 1
      first = .true.
      last  = .false.

C     The SCF loop is executed at least once, 
C     except if nscf <= 0.
      if (nscf <= 0) then
         last = .true.
      endif

#ifdef DEBUG
        call write_debug( '  PRE setup_H0' )
#endif
      call setup_H0()

      DO WHILE (.not. last)
        call timer( 'IterSCF', 1 )
#ifdef DEBUG
        call write_debug( '  PRE setup_hamiltonian' )
#endif

C       Will use H0 computed above
        call setup_hamiltonian( iscf )

#ifdef DEBUG
        call write_debug( '  PRE compute_dm' )
#endif
        call timer( 'compute_dm', 1 )
        call compute_dm( iscf )
        call timer( 'compute_dm', 2 )

#ifdef DEBUG
        call write_debug( '  PRE mixer' )
#endif

        call mixer( iscf, dDmax )

#ifdef DEBUG
        call write_debug( '  PRE scfconvergence_test' )
#endif
        call scfconvergence_test( first, iscf, last, dDmax )

        call timer( 'IterSCF', 2 )
        if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )

        iscf = iscf + 1
        first = .false.

C       The last scheduled iteration is used to compute forces
C       and stresses, hence the exit if iscf==nscf.
C       The ">" clause catches the nscf==1 case. For this,
        if (iscf >= nscf) then
          last = .true.
        endif
      ENDDO

#ifndef RGT
C
C     RGT: Write DM and LWF files
C
#ifdef DEBUG
        call write_debug( '  PRE Write DM' )
#endif
      call iodm( 'write', maxnh, no_l, nspin,
     &                   numh, listhptr, listh, Dscf, found )

#ifdef DEBUG
        call write_debug( '  PRE Write LWF' )
#endif
      if (isolve .eq. 1) then
        call iolwf( 'write', no_s, nbasisCloc, no_l, maxnc, 
     &              found, nspin)
      endif
#endif

C       End of standard SCF loop.
C       Do one more pass to compute forces and stresses
      if (last) then
#ifdef DEBUG
        call write_debug( '  PRE post_scf_work' )
#endif
        call post_scf_work( iscf )
      else
        STOP 'Bad logic in SCF loop'
      endif

#ifdef DEBUG
        call write_debug( '  PRE state_analysis' )
#endif
      call state_analysis( istep )

#if defined(MPI_TRACE) && defined(MPI)
      call mpitrace_eventandcounters( 1000, 0 )
#endif
!--------------------------------------------------------------------------- END
        END subroutine siesta_forces
      END module m_siesta_forces
