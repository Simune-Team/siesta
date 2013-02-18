      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end( )
      use alloc,           only : alloc_report, de_alloc
      use siesta_cml,      only : cml_p, mainXML, cmlEndModule
      use siesta_cmlsubs,  only : siesta_cml_exit
#ifdef DEBUG_XC
      USE siestaXC,        only : closeDebugOutputFile
#endif /* DEBUG_XC */
      use m_timestamp,     only : timestamp
      use m_wallclock,     only : wallclock
      use parallel,        only : IOnode, worker
#ifdef MPI
      use mpi_siesta
#endif
      use alloc,           only : re_alloc
      use m_pulay,         only : resetPulayArrays
      use moreMeshSubs,    only : resetMeshDistr
      use densematrix,     only : resetDenseMatrix
      use sparse_matrices, only : resetSparseMatrices
      use m_fixed,         only : resetFixedPointers
      use band,            only : reset_bands
      use meshphi,         only : resetMeshPhi
      use meshdscf,        only : resetDscfPointers
      use m_dscfcomm,      only : resetdscfComm

      implicit none

#ifdef MPI
      integer MPIerror
#endif
!------------------------------------------------------------------------- BEGIN

      if (worker) then
!     Free allocated memory
      call resetSparseMatrices( )
      call resetPulayArrays( )
      call resetMeshDistr( )
      call resetdscfComm( )
      call resetDenseMatrix( )
      call resetFixedPointers( )
      call reset_bands( )
      call resetMeshPhi( )
      call resetDscfPointers( )

!     Output memory use up to the end of the program
!!      call printmemory( 6, 1 )

!     Print allocation report
      call alloc_report( printNow=.true. )

#ifdef DEBUG_XC
! Close debug outup file
      call closeDebugOutputFile()
#endif /* DEBUG_XC */

!     Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

!     Print final date and time
      if (IOnode) then
        call timestamp('End of run')
        call wallclock('End of run')
      endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif      
      endif  ! worker

!     Finalize MPI
#ifdef MPI
      call MPI_Finalize( MPIerror )
#endif

      if (IOnode) then
         call system("/usr/bin/touch 0_NORMAL_EXIT")
      endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
