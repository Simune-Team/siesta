program getResults

! Collects results of siesta jobs contained in a list, read from standard input
! J.M.Soler. Nov.2012

  use jobList, only: collectResults => getResults

  implicit none
  integer :: unit

  unit = 5
  call collectResults( unit )

end program getResults


