module m_pexsi
  use precision, only: dp
  use iso_c_binding

  integer(c_intptr_t), public :: plan

  public :: pexsi_initialize_scfloop
  public :: pexsi_finalize_scfloop

  private

CONTAINS

subroutine pexsi_initialize_scfloop(World_Comm,npPerPole,mpirank)
  use f_ppexsi_interface
  integer, intent(in) :: npPerPole, mpirank
  integer, intent(in) :: World_Comm

  integer :: numProcRow, numProcCol
  integer :: outputFileIndex, info

numProcRow = sqrt(dble(npPerPole))
numProcCol = numProcRow

if ((numProcRow * numProcCol) /= npPerPole) then
  call die("not perfect square")
endif

outputFileIndex = mpirank

   plan = f_ppexsi_plan_initialize(&
      World_Comm,&
      numProcRow,&
      numProcCol,&
      outputFileIndex,&
      info) 

if (mpirank == 0) then
   print *, "Info in plan_initialize: ", info
endif
end subroutine pexsi_initialize_scfloop

subroutine pexsi_finalize_scfloop() ! mpirank)
  use f_ppexsi_interface
!  integer, intent(in) :: mpirank
  
  integer :: info

  call f_ppexsi_plan_finalize( plan, info )

!  if (mpirank == 0) then
!     print *, "Info in plan_finalize: ", info
!  endif
end subroutine pexsi_finalize_scfloop
  
end module m_pexsi
