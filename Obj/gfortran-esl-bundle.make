#
SIESTA_ARCH=gfortran-esl-bundle-0.4
#
# NOTE: To be used with the ESL bundle V0.4.
#
#  Set the bundle installation directory here
#  Note that it must include libxc
#
BUNDLE_DIR=/Users/ag/lib/Bundle-0.4-with-gxc-0.8.5.1
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
PSML_ROOT=$(BUNDLE_DIR)
XMLF90_ROOT=$(BUNDLE_DIR)
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
 GRIDXC_INCFLAGS:=-I $(BUNDLE_DIR)/include/gridxc_mpi
 GRIDXC_LIBS:=-L $(BUNDLE_DIR)/lib -lgridxc_mpi -lxcf90 -lxc
else
 FC=$(FC_SERIAL)
 GRIDXC_INCFLAGS:=-I $(BUNDLE_DIR)/include/gridxc
 GRIDXC_LIBS:=-L $(BUNDLE_DIR)/lib -lgridxc -lxcf90 -lxc
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
# GRIDXC is treated differently (see above) since the bundle might
# not provide consistent .mk files for it. Note that the bundle must
# include libxc. If not, delete '-lxc90 -lxc' from GRIDXC_LIBS above.
# Note also that the libxc .mod and .h files are found because they
# are in the flat list under $(BUNDLE_DIR)/include, which is searched
# (twice...) since XMLF90_INCFLAGS and PSML_INCFLAGS expand to it.
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
