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
SIESTA_ARCH=sgi64-mpi_fermat
#
# This file seems to work for SGI systems using 64-bit compiled
# Scalapack and Blacs libraries from netlib, and *some version* (perhaps
# SGI's own?) of MPI.
#
# Note that the Scalapack and Blacs library files must be linked from
# their standard places to the building directory... 
#
# Note that the locations of the libraries are configured for the
# Fermat system of the CSAR service at Manchester. Details must
# changed for other machines accordingly.
#
FC=f90 -64
#
FFLAGS=  -O3 
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -O0
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/include
DEFS_MPI=-DMPI
#
LIBS= -lscalapack -lblacs -lpblas -ltools \
      /usr/local/unsupported/numerical/lib64/blacsCinit_MPI-O2K_64-0.a \
      /usr/local/unsupported/numerical/lib64/blacs_MPI-O2K_64-0.a \
      -lscs -lcomplib.sgimath -lmpi 
SYS=bsd
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<
#








