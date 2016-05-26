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
      module files
!
!     Contains the short system label, used to generate file names
!     slabel is currently set in reinit.
!
      integer, parameter, public                  :: label_length = 60
      character(len=label_length), save, public   :: slabel

      ! Derived type to hold some output file names
      type, public:: filesOut_t
        character(len=label_length+6)::
     &    rho     = ' ',  ! (pseudo)electron density
     &    tdrho   = ' ',  ! (pseudo) time-dependent electron density
     &    drho    = ' ',  ! diff. between SCF and atomic electron densities
     &    rhoxc   = ' ',  ! electron density including nonlinear core correction
     &    psch    = ' ',  ! soft diffuse ionic charge
     &    toch    = ' ',  ! total ionic+electronic charge
     &    vh      = ' ',  ! Hartree electrostatic potential
     &    vt      = ' ',  ! total effective potential
     &    vna     = ' '   ! neutral-atom potential
      end type filesOut_t

      private

      end module files
