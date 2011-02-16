      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end( )
      use alloc,           only : alloc_report
      use siesta_cml,      only : cml_p, mainXML, cmlEndModule
      use siesta_cmlsubs,  only : siesta_cml_exit
      use siesta_master,   only : siesta_subroutine  ! Is siesta a subroutine?
#ifdef DEBUG_XC
      USE siestaXC,        only : closeDebugOutputFile
#endif /* DEBUG_XC */
      use m_timestamp,     only : timestamp
      use m_wallclock,     only : wallclock
      use parallel,        only : IOnode
#ifdef MPI
      use mpi_siesta,      only : MPI_Finalize
#endif
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
#ifdef DEBUG
      call write_debug( 'PRE SIESTA_END' )
#endif
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

! Finalize MPI, unless siesta is running as a subroutine
#ifdef DEBUG
      call write_debug( 'POS SIESTA_END' )
#endif
#ifdef MPI
      if (.not.siesta_subroutine) call MPI_Finalize( MPIerror )
#endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif      
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
