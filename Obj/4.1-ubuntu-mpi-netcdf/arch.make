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
# module add mkl/10.1 gcc/4.6.1 scalapack/1.8.0_mpibull2-1.3.9 mpich2/1.4.1-gcc-4.6.1 blacs/1.1-mpibull2-1.3.9
#

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

SIESTA_ARCH=x86_64-mkl_mpi

FPP=
FPP_OUTPUT= 
FC=mpif90 -I. -DMPI  # gfortran v 4.1.2 ???
FC_ASIS=gfortran 
RANLIB=ranlib

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS=-g -O2
FPPFLAGS=-DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DHAVE_MPI -DHAVE_SCALAPACK -DCDF
LDFLAGS= # -static

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

BLAS_LIBS=-lblas
LAPACK_LIBS=-llapack
BLACS_LIBS= # -lblacs-openmpi -lblacsF77init-openmpi
#SCALAPACK_LIBS= -lscalapack-openmpi -lblacs-openmpi -lblacsF77init-openmpi -L${MKLROOT}/lib -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

#SCALAPACK_LIBS= -lscalapack-openmpi -lblacs-openmpi -L${MKLROOT}/lib -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

SCALAPACK_LIBS= -lscalapack-openmpi -lpthread -lm -ldl

NETCDF_ROOT = /usr
NETCDF_INCFLAGS = -I$(NETCDF_ROOT)/include
NETCDF_LIBS = -L$(NETCDF_ROOT)/lib -lnetcdff -lnetcdf
NETCDF_INTERFACE=

LIBS=$(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(NETCDF_LIBS)

#SIESTA needs an F90 interface to MPI
#This will give you SIESTA's own implementation
#If your compiler vendor offers an alternative, you may change
#to it here.
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE= . #/usr/local/software/mpich2-1.4.1-gcc-4.6.1/include

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

