#
SIESTA_ARCH = gfortran-mpi-easybuild
#
#--------------------------------------------------------
# Use these symbols to request particular features
WITH_MPI = 1
WITH_NETCDF = 1
WITH_SPGLIB = 1
#
#--------------------------------------------------------
# Make sure you have the appropriate symbols
# (Either explicitly here, or through shell variables, perhaps
#  set by a module system)
#
NETCDF_ROOT = $(EBROOTNETCDFMINFORTRAN)
PSML_ROOT = $(EBROOTLIBPSML)
XMLF90_ROOT = $(EBROOTXMLF90)
GRIDXC_ROOT = $(EBROOTLIBGRIDXC)
#
#----------------------------------------------------
# Libraries
SCALAPACK_LIBS = -lscalapack
LAPACK_LIBS = -lopenblas
#COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a  # Generic built-in
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC_PARALLEL=mpifort
FC_SERIAL=gfortran
#
FC_ASIS=$(FC)

#FFLAGS = -g -O0 -std=f2003 -fall-intrinsics
FFLAGS = -g -O2 -fbacktrace #-fimplicit-none
FFLAGS_CHECKS = -g -O0 -fcheck=all -Warray-temporaries
#FFLAGS = -O0 -g -fcheck=all  -Warray-temporaries
FFLAGS_DEBUG = -g -O0
RANLIB = :
#
#--------------------------------------------------------
# Nothing should need to be changed below
#
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
ifdef WITH_SPGLIB
 FPPFLAGS_SPG = -DSIESTA__HAS__SPGLIB
endif
#
ifdef WITH_MPI
 FC=$(FC_PARALLEL)
 MPI_INTERFACE=libmpi_f90.a
 MPI_INCLUDE=.      # Note . for no-op
 FPPFLAGS_MPI=-DMPI -DMPI_TIMING
 LIBS +=$(SCALAPACK_LIBS)
else
 FC=$(FC_SERIAL)
endif

LIBS += $(LAPACK_LIBS) $(COMP_LIBS)

SYS = nag
FPPFLAGS = $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) $(FPPFLAGS_SPG) -DF2003 
#
#---------------------------------------------
include $(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk
include $(PSML_ROOT)/share/org.siesta-project/psml.mk
include $(GRIDXC_ROOT)/share/org.siesta-project/gridxc.mk
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
