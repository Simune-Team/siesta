      MODULE MPI__INCLUDE
!
!
!       Michael Hennecke
!       A Fortran 90 interface to MPI version 1.1
!       RZ Uni Karlsruhe, Internal Report 63/96
!
!         Auxiliary module MPI__INCLUDE
!
!         03-Oct-1996, hennecke@rz.uni-karlsruhe.de (v0.9c beta)
!
!       Permission is granted to copy and distribute this file
!       or modified versions of this file for no fee, provided the 
!       copyright notice and this permission notice are preserved 
!       on all copies.
!
!       (C) 1996  Michael Hennecke, RZ Universitaet Karlsruhe
!
!
        IMPLICIT NONE

        INCLUDE 'mpif.h'

        public :: mpi_real
        public :: mpi_complex
        public :: mpi_double_complex
        public :: mpi_double_precision
        public :: mpi_2double_precision
        public :: mpi_integer, mpi_character, mpi_logical
        public :: mpi_maxloc, mpi_sum, mpi_max, mpi_lor
        public :: mpi_status_size
        public :: mpi_comm_world

        private

      END MODULE MPI__INCLUDE
