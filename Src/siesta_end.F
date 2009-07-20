      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end()
      use alloc, only: alloc_report
      use siesta_cml
      USE m_debug,     only : closeDebugOutputFile
      use m_timestamp, only : timestamp
      use m_wallclock, only : wallclock
      use parallel, only: IOnode
#ifdef MPI
      use mpi_siesta, only: MPI_Finalize
#endif
      implicit none

#ifdef MPI
      integer MPIerror
#endif

!------------------------------------------------------------------------- BEGIN
! Output memory use up to the end of the program
      call printmemory( 6, 1 )

! Print allocation report
      call alloc_report( printNow=.true. )

! Close debug outup file
      call closeDebugOutputFile()

! Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

! Print final date and time
      if (IOnode) then
        call timestamp('End of run')
        call wallclock('End of run')
      endif

! Finalize MPI
#ifdef MPI
      call MPI_Finalize( MPIerror )
#endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif      
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
