# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

SIESTA_ARCH=cscs-ibm-blanc-mpi

FPP=/opt/ibmcmp/xlf/10.1/exe/cpp
FPP_OUTPUT= 
FC=mpfort -m64 -compiler xlf95_r
FC_SERIAL=xlf95_r -q64
FFLAGS=-O3 -qstrict -qarch=pwr5 -qtune=pwr5 -qcache=auto -qessl   #-I/apps/netcdf/64/include
# For debugging pourposes
#FC=mpfort -m64 -compiler xlf95_r -qfullpath -qlinedebug -g 
RANLIB=ranlib

SYS=xlf

SP_KIND=1
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FPPFLAGS= -WF,-DMPI -WF,-DFC_HAVE_ABORT -WF,-DCDF
LDFLAGS=

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=-qfixed -qsuffix=cpp=f
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=-qfixed -qsuffix=cpp=F
FPPFLAGS_free_F90=

BLAS_LIBS=-lessl
LAPACK_LIBS=-L/apps/lapack-3/64/lib -llapack_ppc64 -lessl
BLACS_LIBS=/apps/scalapack/64/lib/blacsCinit_MPI-ppc64-0.a /apps/scalapack/64/lib/blacsF77init_MPI-ppc64-0.a /apps/scalapack/64/lib/blacs_MPI-ppc64-0.a
SCALAPACK_LIBS=/apps/scalapack/64/lib/blacsCinit_MPI-ppc64-0.a /apps/scalapack/64/lib/blacs_MPI-ppc64-0.a -L/apps/scalapack/64/lib -lscalapack

COMP_LIBS=dc_lapack.a 

NETCDF_ROOT=/apps/netcdf/64
INCFLAGS=-I$(NETCDF_ROOT)/include
NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf
DEFS_CDF=-DCDF
NETCDF_INTERFACE=

LIBS=$(NETCDF_LIBS) $(SCALAPACK_LIBS) $(BLACS_LIBS) $(LAPACK_LIBS) 

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
