# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
SIESTA_ARCH=osfdxml-mpich
#
# The exact name and location of the libraries might vary.
#
FC=mpif90
FFLAGS= -fast -tune host
FFLAGS_DEBUG= -g3 -fast -tune host
RANLIB=echo
LIBS= /usr/local/BLACS/LIB/blacs_MPI-ALPHA-0.a \
      /usr/local/BLACS/LIB/blacsCinit_MPI-ALPHA-0.a \
      /usr/local/BLACS/LIB/blacsF77init_MPI-ALPHA-0.a \
      /usr/local/SCALAPACK/pblas_ALPHA.a \
      /usr/local/SCALAPACK/scalapack_ALPHA.a \
      /usr/local/BLACS/LIB/blacs_MPI-ALPHA-0.a \
      /usr/local/SCALAPACK/pblas_ALPHA.a \
      /usr/local/SCALAPACK/tools_ALPHA.a \
      -ldxml
SYS=bsd
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
DEFS=-DMPI
MPILIB=
CPP=/bin/cpp -P
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(CPP) $(DEFS) $< > $*.f90
	$(FC) -c $(FFLAGS) $(INCFLAGS) $*.f90
	@rm -f $*.f90
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#

