subroutine memory_all(str,comm)
 use m_rusage, only :  rss_max
#ifdef MPI  
 use mpi
#endif
 character(len=*), intent(in) :: str
 integer, intent(in)          :: comm

 integer :: mpierror
 real    :: max_mem, min_mem, my_mem
 integer :: nprocs, myrank

#ifdef MPI
    call MPI_Comm_Size( Comm, nprocs, MPIerror )
    call MPI_Comm_Rank( Comm, myrank, MPIerror )
#else
    nprocs = 1
    myrank = 0
#endif

 my_mem = rss_max()

#ifdef MPI 
 call MPI_Reduce(my_mem,max_mem,1,MPI_Real,MPI_max,0,comm,MPIerror)
 call MPI_Reduce(my_mem,min_mem,1,MPI_Real,MPI_min,0,comm,MPIerror)
#else
 max_mem = my_mem
 min_mem = my_mem
#endif
 
 if (myrank == 0) then
    write(6,"(a,2f12.2)") " &m -- Peak memory (Mb) " // trim(str) // " (max,min): ",  max_mem, min_mem
 endif

end subroutine memory_all
