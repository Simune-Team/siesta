# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=intel9-mkl8
#
# Intel fortran compiler 9 for linux with mkl 8 optimized blas and lapack
#
# Be sure to experiment with different optimization options.
# You have quite a number of combinations to try...
#
# Note: The -mp1 option is necessary to recover IEEE floating point precision.
#
FC=ifort
#
FFLAGS= -w -xP -O3 -mp1
EXTRA_LIBS=-lpthread -lsvml
FFLAGS_DEBUG= -g 
LDFLAGS= -static
COMP_LIBS=
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#
GUIDE=/opt/intel/mkl/8.0.1/lib/32/libguide.a
LAPACK=/opt/intel/mkl/8.0.1/lib/32/libmkl_lapack.a
BLAS=/opt/intel/mkl/8.0.1/lib/32/libmkl_ia32.a
LIBS=$(LAPACK) $(BLAS)  $(GUIDE) $(EXTRA_LIBS)
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








