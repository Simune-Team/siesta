C****************************************************************
C  Data for passing to function numb                            *
C****************************************************************
      module numbvect

        implicit none

        integer
     .    p, nb

        real*8
     .    qtot

        real*8, dimension(:), allocatable, save ::
     .    c

        real*8, dimension(:,:), allocatable, save ::
     .    rr

      end module numbvect
