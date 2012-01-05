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
# ------------
# Description:
#              Intel compiler/mkl  V10
#              OpenMPI with support for Intel compiler V10
#              MKL V10 libraries, including a version of BLACS provided
#                             by Intel for the openmpi framework,
#                             and Intel's own Scalapack, Lapack, and BLAS.
#
# Execution:
#
#           $(OPENMPI_ROOT)/bin in PATH
#           $(OPENMPI_ROOT)/lib in LD_LIBRARY_PATH
#
SIESTA_ARCH=ompi-intel10-mkl-all
#
#  Edit these parameters for your installation
#
OPENMPI_ROOT=/share/apps/openmpi-1.4.2-intel
MKL_ROOT=/opt/intel/mkl/10.0.3.020/lib/em64t
METIS_LIB=$(HOME)/lib/libmetis.a
#
#--------------------------------------------------------------------------
#
FC=$(OPENMPI_ROOT)/bin/mpif90
#
#  You can play with other optimization options
#  I am not sure whether the compiler attempts to multithread the code
#
# Note: The -mpX option is necessary to recover IEEE floating point precision.
#
FFLAGS= -w  -O2 -mp
FFLAGS_CHECKS= -O0 -g -mp -C -debug all -traceback
FFLAGS_DEBUG= -g 
EXTRA_LIBS=-lpthread -lsvml
COMP_LIBS=
RANLIB=echo
#
NETCDF_ROOT=/share/apps/netcdf-3.6.2-ifort
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
FPPFLAGS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI
#
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
#
INTEL_LIBS_ALL=-L$(MKL_ROOT) \
          -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack \
          -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#
LIBS=$(INTEL_LIBS_ALL)   $(METIS_LIB) $(NETCDF_LIBS)
#
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
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








