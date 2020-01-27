#
SIESTA_ARCH=MN4-ifort
#
# NOTE: To be used with the "last" libgridxc of the 0.8 series,
#       without auto-tools support
#
#--------------------------------------------------------
# Use these symbols to request particular features
# To turn on, set '=1'.
#
WITH_FLOOK=1
WITH_ELSI=1
WITH_EXTERNAL_ELPA_IN_ELSI=
WITH_MPI=1
WITH_NETCDF=1
WITH_NCDF=1
# This will not work until libgridxc 0.9.X
WITH_GRID_SP=
#
#--------------------------------------------------------
# Make sure you have the appropriate symbols
# (Either explicitly here, or through shell variables, perhaps
#  set by a module system)
#
#XMLF90_ROOT=
#PSML_ROOT=
#GRIDXC_ROOT=
#NETCDF_ROOT=/usr/include  
SCALAPACK_LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LAPACK_LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FFTW_ROOT=/apps/FFTW/3.3.6/INTEL/IMPI
#LAPACK_LIBS=-lveclibfort     # Appropriate for MacOS (Homebrew)
#COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a  # Generic built-in
#
# Needed for PEXSI (ELSI) support
LIBS_CPLUS=-lstdc++ 
#--------------------------------------------------------
#
# Define compiler names and flags
#
FC_PARALLEL=mpiifort
FC_SERIAL=ifort
#
FPP = $(FC_SERIAL) -E -P -x c
#
FFLAGS = $(FCFLAGS)
FFLAGS_DEBUG= -g -O0

RANLIB=echo
FC_ASIS=$(FC_SERIAL)
COMP_LIBS=
#
#--------------------------------------------------------
# Nothing should need to be changed below
#
ifdef WITH_GRID_SP
   $(error GRID_SP option does not work with libgridxc < 0.9.X)
endif
#

FPPFLAGS= -DF2003 

LIBS=

ifdef WITH_ELSI
 ifndef ELSI_ROOT
   $(error you need to define ELSI_ROOT in your arch.make)
 endif
 #  Add the second symbol for MAGMA and EigenExa support
 FPPFLAGS_ELSI=-DSIESTA__ELSI # -DSIESTA__ELSI_2_4_SOLVERS

 ELSI_INCFLAGS = -I$(ELSI_ROOT)/include

 ifdef WITH_EXTERNAL_ELPA_IN_ELSI
   ifndef ELPA_ROOT
     $(error you need to define ELPA_ROOT in your arch.make)
   endif
   ifndef ELPA_INCLUDE_DIRECTORY
     # It cannot be generated directly from ELPA_ROOT...
     $(error you need to define ELPA_INCLUDE_DIRECTORY in your arch.make)
   endif
   ELSI_INCFLAGS += -I$(ELPA_INCLUDE_DIRECTORY)
   #ELPA_ROOT=/path/to/external/elpa/installation
 else
   ELPA_ROOT=$(ELSI_ROOT)
 endif
  ELSI_LIB = -L$(ELSI_ROOT)/lib -lelsi \
               -lfortjson -lOMM -lMatrixSwitch \
               -lNTPoly -lpexsi -lsuperlu_dist \
               -lptscotchparmetis -lptscotch -lptscotcherr \
               -lscotchmetis -lscotch -lscotcherr \
               -L$(ELPA_ROOT)/lib -lelpa

  INCFLAGS += $(ELSI_INCFLAGS)
  FPPFLAGS += $(FPPFLAGS_ELSI)
  LIBS +=$(ELSI_LIB) $(LIBS_CPLUS)
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
ifdef WITH_NCDF
 ifndef WITH_NETCDF
   $(error For NCDF you need to define also WITH_NETCDF in your arch.make)
 endif
 FPPFLAGS += -DNCDF -DNCDF_4
 COMP_LIBS += libncdf.a libfdict.a
endif
#
ifdef WITH_FLOOK
 ifndef FLOOK_ROOT
   $(error you need to define FLOOK_ROOT in your arch.make)
 endif
 FLOOK_INCFLAGS=-I$(FLOOK_ROOT)/include
 INCFLAGS += $(FLOOK_INCFLAGS)
 FLOOK_LIBS= -L$(FLOOK_ROOT)/lib -lflookall -ldl
 FPPFLAGS_FLOOK= -DSIESTA__FLOOK
 LIBS +=$(FLOOK_LIBS)
 COMP_LIBS+= libfdict.a
endif
#
ifdef WITH_MPI
 FC=$(FC_PARALLEL)
 MPI_INTERFACE=libmpi_f90.a
 MPI_INCLUDE=.      # Note . for no-op
 FPPFLAGS_MPI=-DMPI -DMPI_TIMING
 LIBS +=$(SCALAPACK_LIBS)
 LIBS +=$(LAPACK_LIBS)
else
 FC=$(FC_SERIAL)
 LIBS += $(LAPACK_LIBS) $(COMP_LIBS)
endif

SYS=nag
FPPFLAGS += $(FPPFLAGS_CDF) $(FPPFLAGS_GRID) $(FPPFLAGS_MPI) $(FPPFLAGS_FLOOK)
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
