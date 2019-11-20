#
SIESTA_ARCH=MN4-psml-elsi-ext-elpa
#
# NOTE: To be used with the "last" libgridxc of the 0.8 series,
#       without auto-tools support
#
#--------------------------------------------------------
# Use these symbols to request particular features
# To turn on, set '=1'.
#
WITH_ELSI=1
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
ELSI_ROOT=/home/pr1eka00/pr1eka01/lib/prace/ifort-17.4/elsi-ext-elpa/2.4.1
ELPA_ROOT=/home/pr1eka00/pr1eka01/lib/prace/ifort-17.4/elpa-2019.05.002
ELPA_INCLUDE_DIRECTORY=/home/pr1eka00/pr1eka01/lib/prace/ifort-17.4/elpa-2019.05.002/include/elpa-2019.05.002/modules
NETCDF_ROOT=/apps/NETCDF/4.4.1.1/INTEL/IMPI
MKLROOT=/apps/INTEL/2017.4/mkl
#
#PSML_ROOT=
#XMLF90_ROOT=
#GRIDXC_ROOT=   (Must be 0.8.5....)
#
LIBS_CPLUS=-lstdc++ 
SCALAPACK_LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LAPACK_LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#
#LAPACK_LIBS=-lveclibfort     # Appropriate for MacOS (Homebrew)
#COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a  # Generic built-in
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC=mpiifort
FC_SERIAL=ifort
#
FC_ASIS=$(FC)

FPP = $(FC) -E -P -x c
#
FFLAGS = -O2 -ip
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
#--------------------------------------------------------
# Nothing should need to be changed below
#
ifdef WITH_GRID_SP
   $(error GRID_SP option does not work with libgridxc < 0.9.X)
endif
#  Add the second symbol for MAGMA and EigenExa support
#
FPPFLAGS_ELSI=-DSIESTA__ELSI # -DSIESTA__ELSI_2_4_SOLVERS
#
ELSI_INCFLAGS = -I$(ELSI_ROOT)/include -I$(ELPA_INCLUDE_DIRECTORY)
ELSI_LIB = -L$(ELSI_ROOT)/lib -lelsi \
               -lfortjson -lOMM -lMatrixSwitch \
               -lNTPoly -lpexsi -lsuperlu_dist \
               -lptscotchparmetis -lptscotch -lptscotcherr \
               -lscotchmetis -lscotch -lscotcherr \
               -L$(ELPA_ROOT)/lib -lelpa

FPPFLAGS= -DF2003 

ifdef WITH_ELSI
 ifndef ELSI_ROOT
   $(error you need to define ELSI_ROOT in your arch.make)
 endif
 ifndef ELPA_INCLUDE_DIRECTORY
   $(error you need to define ELPA_INCLUDE_DIRECTORY in your arch.make)
 endif
 INCFLAGS += $(ELSI_INCFLAGS)
 FPPFLAGS += $(FPPFLAGS_ELSI)
 LIBS +=$(ELSI_LIB)
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
else
 FC=$(FC_SERIAL)
 LIBS += $(LAPACK_LIBS) $(COMP_LIBS)
endif

FC_ASIS=$(FC)

SYS=nag
FPPFLAGS += $(FPPFLAGS_CDF) $(FPPFLAGS_GRID) $(FPPFLAGS_MPI)
#
#---------------------------------------------
include $(XMLF90_ROOT)/share/org.siesta-project/xmlf90.mk
include $(PSML_ROOT)/share/org.siesta-project/psml.mk
include $(GRIDXC_ROOT)/gridxc.mk
#
# We assume that libgridxc
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
