program test_gridxc

use gridxc, only: gridxc_t, cellxc
use gridxc, only: gridxc_initialize

type(gridxc_t) :: gxc

call grixc_initialize(gxc,comm=MPI_Comm_World)
! get density, etc

call setmeshdistr(gxc,....)

call cellxc(gxc,dens,box,.....)

end program test_gridxc

