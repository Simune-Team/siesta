      Program SIESTA

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      USE m_steps, only: inicoor, fincoor
      use alloc

      use parallel, only: SIESTA_worker ! whether part of Siesta's communicator
      use parallel, only: ionode
      use mpi_siesta, only: true_mpi_comm_world
      use m_mpi_utils, only: broadcast

      implicit none

      integer :: istep
      logical :: relaxd

! Notes for PEXSI operation
!
! A subset of nodes carries out non-PEXSI Siesta operations 
! (i.e., setting up H, moving atoms, diagonalizing...). 
! These are tagged as "SIESTA_worker" (admittedly, a bad name)
! All nodes are involved in the PEXSI electronic-structure solver,
! and in the new LocalDOS computation based on selected inversion.
!
! 'siesta_init', 'siesta_forces', and 'siesta_analysis' need to
! be called by all nodes. 
!
! In some cases, the result of a computation by "SIESTA_worker" nodes needs
! to be broadcast to guarantee proper control-flow logic (for example,
! the "relaxd" variable from 'siesta_move'.
!
!----------------------------------------------------------------- BEGIN
!      if (ionode) call memory_snapshot("at start of program")
      call siesta_init()
      if (ionode) call memory_snapshot("after siesta_init")

C     Begin of coordinate relaxation iteration
      relaxd = .false.

      ! Broadcast relevant things for program logic
      ! These were set in siesta_options, called only by "SIESTA_workers".
      call broadcast(inicoor,comm=true_MPI_Comm_World)
      call broadcast(fincoor,comm=true_MPI_Comm_World)

      istep  = inicoor
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
        call siesta_forces( istep )

        if (SIESTA_worker) call siesta_move( istep, relaxd )
        call broadcast(relaxd,comm=true_MPI_Comm_World)

        if (.not. relaxd) then
          istep = istep + 1
        endif
        if (ionode) call memory_snapshot("after geometry step")

      ENDDO

C     End of coordinate-relaxation loop 
      call siesta_analysis( relaxd )
      if (ionode) call memory_snapshot("after siesta_analysis")

      call siesta_end()
      if (ionode) call memory_snapshot("at end of program")

!----------------------------------------------------------------------- END
      END program siesta
