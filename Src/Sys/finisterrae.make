#
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

#
# .make file for Finisterrae at CESGA
# NOTES:
#       The location of the system libraries is implicitly determined
#       by the modules loaded by the user (except netcdf)
#
#
SIESTA_ARCH=intel-mkl-ia64-netcdf-mpi.make

FPP=cp
FPP_OUTPUT=
FC=mpif90.mpich
RANLIB=ranlib

SYS=nag

# These could be made more aggressive, with care
FFLAGS= -w -O2 -mp -cpp 
FFLAGS_DEBUG= -g

ARFLAGS_EXTRA=

GUIDE=-lguide
LAPACK=-lmkl_lapack
BLAS=-lmkl_ipf 
SCALAPACK=-lmkl_scalapack -lmkl_blacs
MKL_LIBS=$(SCALAPACK) $(LAPACK) $(BLAS)  $(GUIDE)  -lpthread 
SYS=nag

NETCDF_ROOT=/opt/cesga/netcdf-3.6.2
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf
FPPFLAGS_CDF=-DCDF

FPPFLAGS= -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT $(FPPFLAGS_CDF)
FPPFLAGS= -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT $(FPPFLAGS_CDF)

METIS_LIBS=$(HOME)/lib/metis/lib/libmetis.a


LIBS=$(MKL_LIBS) $(NETCDF_LIBS) $(METIS_LIBS)

MPI_INCLUDE=.
MPI_INTERFACE=libmpi_f90.a


.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90)  $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

