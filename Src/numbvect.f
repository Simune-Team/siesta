C****************************************************************
C  Data for passing to function numb                            *
C****************************************************************
      module numbvect

        implicit none

        integer
     .    p, nb

        real*8
     .    qtot

        real*8, dimension(:), allocatable ::
     .    c

        real*8, dimension(:,:), allocatable ::
     .    rr

      end module numbvect
