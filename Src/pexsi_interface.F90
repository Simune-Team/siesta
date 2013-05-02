!------------------------
include "pexsi_inertia.h"

   info = 0

end subroutine f_ppexsi_inertiacount_interface

!----------------------
include "pexsi_solve.h"

   integer, parameter :: dp = SELECTED_REAL_KIND(10,100)

   integer :: ierr, mpirank

   call MPI_Comm_Rank(comm_global, mpirank, ierr)

   if (mpirank == 0) then
      write(*,*) "PEXSI:"
      write(*,*) "nrows: ", nrows
      write(*,*) "nnz: ", nnz
      write(*,*) "nnzLocal (0): ", nnzLocal
      write(*,*) "numcolLocal (0): ", numColLocal
      write(*,*) "numPole: ", numPole
      write(*,*) "temperature (Ry): ", temperature
      write(*,*) "numElectronExact: ", numElectronExact
      write(*,*) "numElectronTolerance: ", PEXSInumElectronTolerance
      write(*,*) "gap: ", gap
      write(*,*) "DeltaE: ", deltaE
      write(*,*) "muMin: ", muMinInertia
      write(*,*) "muMax: ", muMaxInertia
      write(*,*) "muMaxIter: ", muMaxIter
      write(*,*) "npPerPole: ", npPerPole
      write(*,*) "PoleTolerance: ", poleTolerance
      write(*,*) "ordering: ", ordering
   endif

!  Sample dummy outputs

   info = 0

   muIter = min(muMaxIter,2)

   numElectron = 24.0_dp
   mu = -0.3_dp

   numElectronList(1:muIter-1) = 23.9_dp
   numElectronList(muIter) = numElectron

   muList(1:muIter-1) = -0.29_dp
   muList(muIter) = mu

   numElectronDrvList(1:muIter) = 0.04_dp

   muZeroT = -0.28_dp

end subroutine f_ppexsi_solve_interface

