program countJobs

! Counts siesta jobs contained in a list, read from standard input
! J.M.Soler. Nov.2012

  use jobList, only: count => countJobs

  implicit none
  integer :: nJobs, nLists, unit

  unit = 5
  call count( unit, nJobs, nLists )

  print*,'number of jobs and lists =', nJobs, nLists

end program countJobs


