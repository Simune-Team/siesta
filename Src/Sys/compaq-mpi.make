# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile include file for parallel Compaq machine at Grenoble

SIESTA_ARCH=compaq-mpi

FC=f90 

FFLAGS=  -O2
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -Rabc -ei

NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=

MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
FPPFLAGS_MPI=-DMPI

LIBS = -L/usr/local/lib/scalapack -lscalapack -lpblas -ltools -lblacsF77 \
        -lblacs -lblacsF77 \
        -lblacs -ldxml -lfmpi -lmpi -lelan

SYS=bsd
RANLIB=echo
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) 

# Actual compilation recipes for siesta code.

.F.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)   $(FPPFLAGS) $<
.f90.o:
        $(FC) -c $(FFLAGS) $(INCFLAGS)   $<

