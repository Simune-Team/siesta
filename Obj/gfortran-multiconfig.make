#
SIESTA_ARCH=gfortran-multiconfig
#
# NOTE: To be used with libgridxc v 0.9.3 and above
# when configured with the 'multiconfig' option.
#--------------------------------------------------------
# Use these symbols to request particular features
# To turn on, set '=1'.
#
WITH_MPI=1
WITH_NETCDF=1
WITH_GRID_SP=
#
#--------------------------------------------------------
# Make sure you have the appropriate symbols
# (Either explicitly here, or through shell variables, perhaps
#  set by a module system)
#
#NETCDF_ROOT=/usr/local
#PSML_ROOT=$(HOME)/lib/gfortran-5.2.0/libpsml-1.1.7
#XMLF90_ROOT=$(HOME)/lib/gfortran-5.2.0/xmlf90-1.5.4
#GRIDXC_ROOT=$(HOME)/lib/gfortran-5.2.0/gridxc-0.8.5
#
#----------------------------------------------------
#SCALAPACK_LIBS=....
#LAPACK_LIBS=-lveclibfort     # Appropriate for MacOS (Homebrew)
#COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a  # Generic built-in
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC_PARALLEL=mpif90
FC_SERIAL=gfortran
#

#FFLAGS= -g -O0 -std=f2003 -fall-intrinsics
FFLAGS= -O2 -fbacktrace #-fimplicit-none
#FFLAGS= -O0 -g -fcheck=all #-Warray-temporaries
FFLAGS_CHECKS= -O0 -g -fcheck=all -Warray-temporaries
#FFLAGS= -O0 -g -fcheck=all  -Warray-temporaries
FFLAGS_DEBUG= -g -O0
RANLIB=echo
#
#--------------------------------------------------------
# Nothing should need to be changed below
#
ifdef WITH_GRID_SP
  GRIDXC_CONFIG_PREFIX=sp
  FPPFLAGS_GRID= -DGRID_SP
else
  GRIDXC_CONFIG_PREFIX=dp
endif

ifdef WITH_NETCDF
 ifndef NETCDF_ROOT
   $(error you need to define NETCDF_ROOT in your arch.make)
 endif
 NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
 NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdff
 FPPFLAGS_CDF= -DCDF
 LIBS +=$(NETCDF_LIBS)
endif
#
ifdef WITH_MPI
 FC=$(FC_PARALLEL)
 MPI_INTERFACE=libmpi_f90.a
 MPI_INCLUDE=.      # Note . for no-op
 FPPFLAGS_MPI=-DMPI -DMPI_TIMING
 LIBS +=$(SCALAPACK_LIBS)
 GRIDXC_CONFIG_PREFIX:=$(GRIDXC_CONFIG_PREFIX)_mpi
else
 FC=$(FC_SERIAL)
endif

LIBS += $(LAPACK_LIBS) $(COMP_LIBS)
FC_ASIS=$(FC)

SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_GRID) $(FPPFLAGS_MPI)  -DF2003 
#
#---------------------------------------------
include $(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk
include $(PSML_ROOT)/share/org.siesta-project/psml.mk
include $(GRIDXC_ROOT)/share/org.siesta-project/gridxc_$(GRIDXC_CONFIG_PREFIX).mk
#---------------------------------------------
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
