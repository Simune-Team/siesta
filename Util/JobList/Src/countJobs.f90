program countJobs

! Counts siesta jobs contained in a list, read from standard input
! J.M.Soler. Nov.2012

  use jobList, only: count => countJobs

  implicit none
  integer :: nCores, nJobs, nLists, unit

  unit = 5
  call count( unit, nLists, nJobs, nCores )

  print*,'number of lists, jobs, and cores =', nLists, nJobs, nCores

end program countJobs


