#
WITH_MPI=1
WITH_NETCDF=1
#
SIESTA_ARCH=gfortran-macosx64-openmpi
# The only thing you should change is the location of the libraries
# on your computer
#
FC_PARALLEL=mpif90
FC_SERIAL=gfortran
#
FC_ASIS=$(FC)

#FFLAGS= -g -O0 -std=f2003 -fall-intrinsics
FFLAGS= -O2 -fbacktrace #-fimplicit-none
FFLAGS_CHECKS= -O0 -g -fcheck=all -Warray-temporaries
#FFLAGS= -O0 -g -fcheck=all  -Warray-temporaries
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
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
ifdef WITH_MPI
 FC=$(FC_PARALLEL)
 MPI_INTERFACE=libmpi_f90.a
 MPI_INCLUDE=.      # Note . for no-op
 FPPFLAGS_MPI=-DMPI -DMPI_TIMING
 LIBS +=$(SCALAPACK_LIBS)
else
 FC=$(FC_SERIAL)
endif

LIBS += -lveclibfort

SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)  -DF2003 
#
#---------------------------------------------
include $(XMLF90_ROOT)/share/siesta/xmlf90.mk
include $(PSML_ROOT)/psml.mk
include $(GRIDXC_ROOT)/gridxc.mk
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
