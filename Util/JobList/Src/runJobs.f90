program runJobs

! Runs siesta jobs contained in a list, read from standard input
! J.M.Soler. Nov.2012

  use jobList, only: runFile => runJobs

  implicit none
  integer :: unit

  unit = 5
  call runFile( unit )

end program runJobs


