! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end( )
      use fdf,             only : fdf_shutdown
      use alloc,           only : alloc_report
      use siesta_cml,      only : cml_p, mainXML, cmlEndModule
      use siesta_cmlsubs,  only : siesta_cml_exit
      use siesta_master,   only : siesta_subroutine  ! Is siesta a subroutine?
#ifdef DEBUG_XC
      USE siestaXC,        only : closeDebugOutputFile
#endif /* DEBUG_XC */
      use m_timestamp,     only : timestamp
      use m_wallclock,     only : wallclock, wallclock_shutdown
      use parallel,        only : IOnode, Node
#ifdef MPI
      use mpi_siesta,      only : MPI_Finalize, MPI_SUCCESS
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
      use m_rhog,          only : resetRhog
      use m_io,            only : io_assign, io_close

      implicit none

      integer :: iu
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
      call resetRhog()

!     Output memory use up to the end of the program
!!      call printmemory( 6, 1 )

!     Print allocation report
      call alloc_report( printNow=.true. , shutdown=.true.)

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
        ! Ensure clock file is closed
        call wallclock_shutdown()
      endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif      

! Finalize MPI, unless siesta is running as a subroutine
#ifdef DEBUG
      call write_debug( 'POS SIESTA_END' )
#endif
      !
      ! Add stamp to MESSAGES file
      call message('INFO','Job completed')

      ! Clean the fdf-files
      call fdf_shutdown()

#ifdef MPI
      MPIerror = MPI_SUCCESS
      if (.not.siesta_subroutine) call MPI_Finalize( MPIerror )
#endif

      if (IOnode) then
         call io_assign(iu)
         open(iu, file='0_NORMAL_EXIT', form='FORMATTED')
#ifdef MPI
         write(iu,'(a,i0)') 'SIESTA completed, MPI exit: ',
     &        MPIerror
#else
         write(iu,'(a)') 'SIESTA completed successfully'
#endif
         call io_close(iu)
      endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
