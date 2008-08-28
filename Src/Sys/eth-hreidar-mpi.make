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
 
SIESTA_ARCH=hreidar--Portland-mpich
FPP=
FPP_OUTPUT=
FC=/usr/local/apli/mpich/bin/mpif90
FC_ASIS=$(FC)
MPIFC=/usr/local/apli/mpich/bin/mpif90
RANLIB=echo

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS=-fastsse
FFLAGS_DEBUG= -g -O0
FPPFLAGS= -DMPI
LDFLAGS=

ARFLAGS_EXTRA=

MPI_ROOT_MPICH=/usr/local/apli/mpich
MPI_LIBS_MPICH= -L$(MPI_ROOT_MPICH)/lib -lmpich
#
MPI_ROOT=$(MPI_ROOT_MPICH)
MPI_LIBS=$(MPI_LIBS_MPICH)
LIBS= -L$(HOME)/scalapack-1.8.0  -lscalapack \
       -L$(HOME)/BLACS/LIB -lblacs -lblacsF77init -lblacsCinit -lblacsF77init -lblacs \
       -L/usr/local/apli/pgi5.2/linux86-64/5.2/lib/ -llapack -lblas \
       $(MPI_LIBS)  $(NETCDF_LIBS)
 
COMP_LIBS=dc_lapack.a
#
SYS=cpu_time
DEFS= $(DEFS_CDF) $(DEFS_MPI)
NETCDF_LIBS=
NETCDF_INTERFACE=

#SIESTA needs an F90 interface to MPI
#This will give you SIESTA's own implementation
#If your compiler vendor offers an alternative, you may change
#to it here.
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/apli/mpich/include
DEFS_MPI=-DMPI
DEFS= $(DEFS_CDF) $(DEFS_MPI)

# Compile atom.f without optimization.
#
atom.o:
	$(FC) -c $(FFLAGS_DEBUG) atom.f
#
electrostatic.o:
	$(FC) -c $(FFLAGS_DEBUG) electrostatic.f

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
