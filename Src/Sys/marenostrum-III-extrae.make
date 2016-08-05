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
#              Intel compiler/mkl  on MN
#              OpenMPI with support for Intel compiler
#              MKL libraries, including a version of BLACS provided
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
PEXSI_VERSION=release_mn_v0.5.4
SIESTA_ARCH=intel-mn-openmpi-pexsi-$(PEXSI_VERSION)
#
#
FC=mpif90
#
#  You can play with other optimization options
#  I am not sure whether the compiler attempts to multithread the code
#
#FFLAGS=-g -O0 -debug full -traceback -C
FFLAGS= -g -w -O2 
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
#NETCDF_ROOT=/apps/NETCDF/4.1.3
#NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
#FPPFLAGS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI
#
#NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
#METIS_LIB=$HOME/lib/metis-4.0/libmetis.a
#
# From the "Intel advisor"
#
MKLPATH=/gpfs/apps/MN3/INTEL/mkl/lib/intel64
SUGGESTED_LIBS=$(MKLPATH)/libmkl_scalapack_lp64.a \
               -Wl,--start-group \
                  $(MKLPATH)/libmkl_intel_lp64.a \
                  $(MKLPATH)/libmkl_sequential.a \
                  $(MKLPATH)/libmkl_core.a \
                  $(MKLPATH)/libmkl_blacs_openmpi_lp64.a \
               -Wl,--end-group \
               -lpthread 
#
# Extended interface
#PEXSI_DIR=$(HOME)/code/pexsi
LIN=$(HOME)/lib
#
PEXSI_DIR=$(LIN)
DSUPERLU_DIR=$(LIN)
PARMETIS_DIR=$(LIN)

METIS_LIB        = ${PARMETIS_DIR}/libmetis.a
PARMETISLIB	     = ${PARMETIS_DIR}/libparmetis.a
#SELINV_LIB       = ${PEXSI_DIR}/libselinv.a
DSUPERLU_LIB     = ${DSUPERLU_DIR}/libsuperlu_dist_3.2.a
PEXSI_LIB        = ${PEXSI_DIR}/libpexsi_$(PEXSI_VERSION).a
LIN_LIBS         = ${PEXSI_LIB} ${DSUPERLU_LIB} ${SELINV_LIB} ${PARMETISLIB} ${METIS_LIB} 
#
EXTRAE_LIBS=/apps/CEPBATOOLS/extrae/latest/openmpi/64/lib/libmpitrace.so
#
LIBS=$(LIN_LIBS) $(SUGGESTED_LIBS) $(NETCDF_LIBS) $(EXTRAE_LIBS) -lstdc++ # $(METIS_LIB) 
#
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) -DTRACING
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








