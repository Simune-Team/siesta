program runJobs

! Runs siesta jobs contained in a list
! J.M.Soler. Nov.2012

  use jobList, only: runFile

  implicit none
  integer :: nJobs, nLists, unit

! Read and run the job list from standard input
  unit = 5
  call runFile( unit, nJobs, nLists )
  print'(a,2i6)', 'submitted jobs, lists =', nJobs, nLists

end program runJobs


