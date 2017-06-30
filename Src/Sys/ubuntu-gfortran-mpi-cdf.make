# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#

#  For Ubuntu, with the libraries set up as described by Mark Calleja in
#    http://pelios.csx.cam.ac.uk/~mc321/siesta.html
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

SIESTA_ARCH=x86_64-unknown-linux-gnu--Gfortran

FPP=
FPP_OUTPUT= 
FC=mpif90
RANLIB=echo

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

# Add any other sensible compilation flags here
FFLAGS=-g -O2

FPPFLAGS= -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DCDF -DGRID_DP -DPHI_GRID_SP
LDFLAGS=

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

INCFLAGS=-I/usr/include -I. -I/usr/lib/openmpi/include
BLAS_LIBS=-lblas
LAPACK_LIBS=/usr/lib/lapack/liblapack.a
BLACS_LIBS=-lblacsF77init-openmpi -lblacsCinit-openmpi -lblacs-openmpi
SCALAPACK_LIBS=-lscalapack-openmpi

COMP_LIBS=dc_lapack.a 

NETCDF_LIBS=-lnetcdff -lnetcdf
NETCDF_INTERFACE=libnetcdf_f90.a

LIBS= $(NETCDF_LIBS) -lpthread $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

#SIESTA needs an F90 interface to MPI
#This will give you SIESTA's own implementation
#If your compiler vendor offers an alternative, you may change
#to it here.
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.

#Dependency rules are created by autoconf according to whether
#discrete preprocessing is necessary or not.
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

