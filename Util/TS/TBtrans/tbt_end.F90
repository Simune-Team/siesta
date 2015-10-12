subroutine tbt_end()

#ifdef MPI
  use mpi_siesta, only : MPI_Finalize
#endif
  use parallel, only : Node
  use alloc, only   : alloc_report
  use m_timestamp, only : timestamp
  use m_wallclock, only : wallclock

#ifdef MPI
  integer :: MPIerror
#endif

  ! Stop time counter
  call timer( 'tbtrans', 2 )
  call timer( 'all', 3 )

  ! Print allocation report
  call alloc_report( printNow=.true. )

  if ( Node == 0 ) then
     call timestamp("End of run")
     call wallclock("End of run")
  end if

#ifdef MPI
  call MPI_Finalize(MPIerror)
#endif

end subroutine tbt_end
