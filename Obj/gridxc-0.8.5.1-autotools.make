#
SIESTA_ARCH=gfortran-gridxc-0.8.5.1
#
# NOTE: To be used with the "last" libgridxc of the 0.8 series,
# which is scheduled to be released with the 0.4 ESL bundle after
# being fitted with an auto-tools building system. It has a
# version id of 0.8.5.1.
#
# It still has a not-optimal tree organization (as compared
# with the forthcoming 0.9 series)
#
#--------------------------------------------------------
# Use these symbols to request particular features
# To turn on, set '=1'.
#
WITH_MPI=1
WITH_NETCDF=1
# This will not work until libgridxc 0.9.X
WITH_GRID_SP=
#
#--------------------------------------------------------
# Make sure you have the appropriate symbols
# (Either explicitly here, or through shell variables, perhaps
#  set by a module system)
#
#NETCDF_ROOT=/usr/local
#PSML_ROOT=
#XMLF90_ROOT=
#GRIDXC_ROOT=   (Must be 0.8.5.1...)
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
   $(error GRID_SP option does not work with libgridxc < 0.9.X)
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
 GRIDXC_INCFLAGS:=-I $(GRIDXC_ROOT)/include/gridxc_mpi
 GRIDXC_LIBS:=-L $(GRIDXC_ROOT)/lib -lgridxc_mpi -lxcf90 -lxc
else
 FC=$(FC_SERIAL)
 GRIDXC_INCFLAGS:=-I $(GRIDXC_ROOT)/include/gridxc
 GRIDXC_LIBS:=-L $(GRIDXC_ROOT)/lib -lgridxc -lxcf90 -lxc
endif

LIBS += $(LAPACK_LIBS) $(COMP_LIBS)
FC_ASIS=$(FC)

SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_GRID) $(FPPFLAGS_MPI)  -DF2003 
#
#---------------------------------------------
include $(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk
include $(PSML_ROOT)/share/org.siesta-project/psml.mk
#
# GRIDXC is treated differently (see above) since we do not have
# fully consistent .mk files for it in 0.8.5.1. We assume that libgridxc
# includes libxc. If not, delete '-lxc90 -lxc' from GRIDXC_LIBS above.
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
