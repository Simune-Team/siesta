MODULE Kpoint_grid
  USE precision, only : dp
  implicit none

  integer  :: nkpnt
  integer  :: maxk   = 0
  real(dp) :: kcutof = 0.0_dp


  integer,  dimension(3,3) :: kscell = 0
  real(dp), dimension(3)   :: kdispl = 0.0_dp

  real(dp), allocatable :: kweight(:)
  real(dp), allocatable :: kpoint(:,:)

  CONTAINS

  subroutine setup_Kpoint_grid( ucell )
    implicit none
    real(dp) :: ucell(3,3)

    ! Estimate the actual value of nkpnt
    call kgridinit( ucell, kscell, kdispl, kcutof, nkpnt )

    IF (maxk .eq. 0) then
      maxk = nkpnt
      allocate(kpoint(3,maxk))
      call memory('A','D',3*maxk,'siesta')
      allocate(kweight(maxk))
      call memory('A','D',maxk,'siesta')

      kpoint(1:3,1) = 0.0_dp
      kweight(1)    = 1.0_dp

    ELSE
      ! If number of k points is greater than the previous one - re-size arrays
      if (nkpnt .gt. maxk) then
        maxk = nkpnt

        call memory('D','D',size(kpoint),'siesta')
        deallocate(kpoint)
        call memory('D','D',size(kweight),'siesta')
        deallocate(kweight)

        allocate(kpoint(3,maxk))
        call memory('A','D',3*maxk,'siesta')
        allocate(kweight(maxk))
        call memory('A','D',maxk,'siesta')

        kpoint(1:3,1) = 0.0_dp
        kweight(1)    = 1.0_dp
       endif
    ENDIF

    call kgrid( ucell, kscell, kdispl, nkpnt, kpoint, kweight )

  end subroutine setup_Kpoint_grid

END MODULE Kpoint_grid
