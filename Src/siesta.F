      Program SIESTA

      use m_siesta_init
      use m_siesta_analysis
      use m_siesta_move
      use m_siesta_end
      use m_siesta_forces
      USE m_steps, only: inicoor, fincoor
      use mpi_siesta
      use alloc
      use parallel, only: worker

      implicit none

      integer :: istep, mpierr
      logical :: relaxd

!----------------------------------------------------------------- BEGIN
      call siesta_init()
!      if (.not. worker) goto 999

C     Begin of coordinate relaxation iteration
      relaxd = .false.
      ! Broadcast relevant things for program logic
      call MPI_Bcast(inicoor,1,MPI_integer,0,
     $               true_MPI_Comm_World,MPIerr)
      call MPI_Bcast(fincoor,1,MPI_integer,0,
     $               true_MPI_Comm_World,MPIerr)
      istep  = inicoor
      DO WHILE ((istep.le.fincoor) .AND. (.not. relaxd))
        call siesta_forces( istep )

        if (worker) call siesta_move( istep, relaxd )
      !
        call MPI_Bcast(relaxd,1,MPI_logical,0,
     $                 true_MPI_Comm_World,MPIerr)

        if (.not. relaxd) then
          istep = istep + 1
        endif

      ENDDO

C     End of coordinate-relaxation loop 
      if (worker) call siesta_analysis( relaxd )

      call siesta_end()
!----------------------------------------------------------------------- END
 999  continue
      END program siesta
