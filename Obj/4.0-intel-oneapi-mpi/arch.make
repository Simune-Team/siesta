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
# sudo apt install libmkl-dev intel-mkl-doc libomp-7-doc tbb-examples libtbb-doc
# sudo apt install libmkl-full-dev
# sudo apt install libmkl-scalapack-lp64 
# sudo apt install libmkl-blacs-intelmpi-lp64
# sudo apt install intel-oneapi-mpi-devel
# sudo apt install intel-oneapi-mkl-common-2021.2.0
# sudo apt install intel-oneapi-mkl-devel
# sudo apt install intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic
# sudo apt install intel-oneapi-compiler-fortran
#
# source /opt/intel/oneapi/setvars.sh 
#

.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

SIESTA_ARCH=x86_64-mkl_mpi

FPP=
FPP_OUTPUT= 
FC=mpiifort -I.
FC_ASIS=ifort
RANLIB=ranlib

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS=-g -O2 -mkl=sequential 
FPPFLAGS=-DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DMPI
LDFLAGS= # -static

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

BLAS_LIBS=
LAPACK_LIBS=
BLACS_LIBS= # -lblacs-openmpi -lblacsF77init-openmpi

#SCALAPACK_LIBS= -lscalapack-openmpi -lblacs-openmpi -lblacsF77init-openmpi -L${MKLROOT}/lib -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

#SCALAPACK_LIBS= -lscalapack-openmpi -lblacs-openmpi -L${MKLROOT}/lib -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

SCALAPACK_LIBS=-mkl -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl 


NETCDF_LIBS=
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

