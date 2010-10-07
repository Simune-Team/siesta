MODULE timer_mpi_m

! Disconnectable interface to siesta's timer routine. 
! J.M.Soler. May.2009

CONTAINS

  SUBROUTINE timer_mpi( name, opt )
    character(len=*), intent(in):: name
    integer,          intent(in):: opt

#ifdef MPI_TIMING
    external timer
    call timer( name, opt )
#endif

  END SUBROUTINE timer_mpi

END MODULE timer_mpi_m
