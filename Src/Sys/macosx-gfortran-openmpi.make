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
SIESTA_ARCH=macosx-openmpi
#
#
FC=/opt/openmpi.gfortran/bin/mpif90
#
FC_ASIS=$(FC)
#
#FFLAGS=-O2
FFLAGS= -g -Wall -Wextra -O0 -g -fbacktrace -fbounds-check -Wl,-commons,use_dylibs
FFLAGS_DEBUG= -g -Wall -Wextra -O0 -g -fbacktrace -fbounds-check
LDFLAGS=
RANLIB=echo
LIBS=  
SYS=nag
#
# --- Edit the location of your netcdf files
#
#NETCDF_ROOT=$(HOME)/lib/netcdf-3.6.2-gfortran
#INCFLAGS=-I$(NETCDF_ROOT)/include
#DEFS_CDF=-DCDF
#
DEFS=-DGFORTRAN  -DFC_HAVE_FLUSH -DFC_HAVE_ABORT    # Note this !!
COMP_LIBS=
#
#NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
DEFS_MPI=-DMPI
#
# Specify full paths to avoid picking up the veclib versions..
#
LAPACK=-framework veclib
BLAS=-framework veclib
BLACS=-L$(HOME)/lib/gfortran/openmpi-1.4.2 -lblacsF77 -lblacsC -lblacs
SCALAPACK=-L$(HOME)/lib/gfortran -lscalapack $(BLACS)
METIS=/Users/ag/lib/metis/lib64/libmetis.a

DEFS:=$(DEFS_MPI) $(DEFS_CDF) $(DEFS)
LIBS= $(SCALAPACK) $(LAPACK) $(BLAS) $(NETCDF_LIBS) $(METIS)


DEFS:=$(DEFS_MPI) $(DEFS_CDF) $(DEFS)
LIBS= $(SCALAPACK) $(LAPACK) $(BLAS) $(NETCDF_LIBS) $(METIS)

#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
