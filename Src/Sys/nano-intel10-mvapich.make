# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=nano-intel-mvapich
#
# To run in parallel with InfiniBand, make sure that
# /usr/mpi/intel/mvapich-1.0.1/bin is in your path, and use
# the following incantation to run (in multiples of 8 procs)
#
#  PATH=/usr/mpi/intel/mvapich-1.0.1/bin
#  mpirun_rsh -np $NSLOTS -hostfile $PBS_NODEFILE  siesta < FILE.fdf > OUT
#
#
#--------------------------------------------------------------------------
# Note: The -mpX option is necessary to recover IEEE floating point precision.
#
FC=/usr/mpi/intel/mvapich-1.0.1/bin/mpif90
#
#  You can play with other optimization options
#  I am not sure whether the compiler attempts to multithread the code
#
FFLAGS= -w -O2 -mp
##FFLAGS= -w  -O0 -g -debug full -traceback -C
EXTRA_LIBS=-lpthread -lsvml
FFLAGS_DEBUG= -g 
LDFLAGS=
COMP_LIBS=
RANLIB=echo
#
NETCDF_ROOT=/share/apps/netcdf-3.6.2-ifort
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
FPPFLAGS_CDF=-DCDF
DEFS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
DEFS_MPI=-DMPI
#
METIS_LIBS=/share/apps/metis-4.0/libmetis.a
#
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
MKLPATH=/opt/intel/mkl/10.0.3.020/lib/em64t
BLAS_LIBS=$(MKLPATH)/libmkl_em64t.a
LAPACK_LIBS=$(MKLPATH)/libmkl_lapack.a
BLACS_LIBS=$(MKLPATH)/libmkl_blacs_lp64.a
SCALAPACK_LIBS=$(MKLPATH)/libmkl_scalapack_lp64.a

LINO_LIBS=-L/opt/intel/mkl/10.0.3.020/lib/em64t \
          -lmkl_scalapack_lp64 -lmkl_blacs_lp64 -lmkl_lapack \
          $(METIS_LIBS) $(NETCDF_LIBS) \
          -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

FDN_LIBS=$(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) \
         $(METIS_LIBS) $(NETCDF_LIBS) -liomp5 -pthread
#
LIBS=$(LINO_LIBS)
SYS=nag
DEFS= $(DEFS_CDF) $(DEFS_MPI) ### -DDEBUG
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








