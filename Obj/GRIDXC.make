# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
SIESTA_ARCH=gfortran-macosx64-openmpi
# The only thing you should change is the location of the libraries
# on your computer
#
FC=mpif90
FC_SERIAL=gfortran
#
GRIDXC_ROOT=$(HOME)/lib/gfortran/gridxc-mpi
GRIDXC_SERIAL_ROOT=$(HOME)/lib/gfortran/gridxc
#
FC_ASIS=$(FC)
#
###FFLAGS= -g -O2 # -Wall
FFLAGS_CHECKS= -O0 -g -fcheck=all
FFLAGS= -O0 -g -fcheck=all
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
NETCDF_ROOT=/usr/local
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
FPPFLAGS_CDF=-DCDF
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdff
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI

LIBS=-L/opt/scalapack/openmpi-1.6.1-gfortran/lib \
        -lscalapack -ltmg -lreflapack -lrefblas \
      $(NETCDF_LIBS)

SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) -DF2003
#
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
