! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2006.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
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
