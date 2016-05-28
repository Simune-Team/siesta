# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# ------------
# Description:
#              Intel compiler/mkl  V11 on atto
#              OpenMPI with support for Intel compiler V11
#              MKL V11 libraries, including a version of BLACS provided
#                             by Intel for the openmpi framework,
#                             and Intel's own Scalapack, Lapack, and BLAS.
#
# Execution:
#
#           $(OPENMPI_ROOT)/bin in PATH
#           $(OPENMPI_ROOT)/lib in LD_LIBRARY_PATH
#
#           mpirun -np NPROCS siesta ....
#
SIESTA_ARCH=intel11-openmpi
#
#
FC=/share/apps/openmpi-1.4.3-intel/bin/mpif90
#
#  You can play with other optimization options
#  I am not sure whether the compiler attempts to multithread the code
#
FFLAGS= -w  -O2 -mp
FFLAGS_CHECKS=-g -O0 -debug full -traceback -C
FFLAGS_DEBUG= -g 
LDFLAGS= 
COMP_LIBS=
RANLIB=echo
#
# You might want to turn off FoX for Intel11
#
DUMMY_FOX=--enable-dummy
#
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
METIS_LIB=/share/apps/metis-4.0/libmetis.a
#
# From the "Intel advisor"
#
MKLPATH=/opt/intel/Compiler/11.1/064/mkl/lib/em64t
SUGGESTED_LIBS=$(MKLPATH)/libmkl_scalapack_lp64.a \
               -Wl,--start-group \
                  $(MKLPATH)/libmkl_intel_lp64.a \
                  $(MKLPATH)/libmkl_sequential.a \
                  $(MKLPATH)/libmkl_core.a \
                  $(MKLPATH)/libmkl_blacs_openmpi_lp64.a \
               -Wl,--end-group \
               -lpthread 
#
LIBS=$(SUGGESTED_LIBS) $(NETCDF_LIBS) $(METIS_LIB)
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








