! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module psop_params

      implicit none 
!
!    Hard-wired parameters to be eliminated in the future
!

C INTEGER  NKBMX    : Maximum number of Kleinman-Bylander projectors
C                     for each angular momentum

         integer, parameter, public  :: nkbmx  =    2  

C INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
C                      projectors.

         integer, parameter, public  :: lmaxd  =    4  

C INTEGER  NRMAX    : Maximum number of points in the functions read
C                     from file '.vps' or '.psf' (this number is
C                     determined by the parameter nrmax in the
C                     program atm, which generates the files with
C                     the pseudopotential information). The number
C                     of points in the grid can be redefined if the
C                     pseudopotential is reparametrized.
C                     nrmax = 20000 is a typical safe value in this case.
C                     

         integer, parameter, public  :: nrmax  = 20000

      end module psop_params
