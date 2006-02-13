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
SIESTA_ARCH=altix-intel
#
# Intel fortran compiler for linux with mkl optimized blas and lapack
#
# Be sure to experiment with different optimization options.
# You have quite a number of combinations to try...
#
FC=ifort -ftz
#
FFLAGS= -w -mp1 -O3 -fpe3
FFLAGS_DEBUG= -g -mp -O0 -CB
LDFLAGS=-Vaxlib 
COMP_LIBS=
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
LIBS=$(MKL_LIBS) -lsdsm -lscs -lmpi -ldl
MKL_HOME=/opt/intel-mkl/7.2.1.003/mkl721/lib/64/
MKL_LIBS=-L/opt/mpt/1.12-sgi402r1/lib
SYS=bsd
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#
