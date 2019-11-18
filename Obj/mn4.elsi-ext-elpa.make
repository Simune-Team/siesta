#
SIESTA_ARCH=ifort-MN
# The only thing you should change is the location of the libraries
# on your computer
#
ELSI_ROOT=/home/pr1eka00/pr1eka01/lib/prace/ifort-17.4/elsi-ext-elpa/2.4.1
ELPA_ROOT=/home/pr1eka00/pr1eka01/lib/prace/ifort-17.4/elpa-2019.05.002
ELPA_INCLUDE_DIRECTORY=/home/pr1eka00/pr1eka01/lib/prace/ifort-17.4/elpa-2019.05.002/include/elpa-2019.05.002/modules
NETCDF_ROOT=/apps/NETCDF/4.4.1.1/INTEL/IMPI
MKLROOT=/apps/INTEL/2017.4/mkl
#
FC=mpiifort
FC_SERIAL=ifort
#
FC_ASIS=$(FC)

FPP = $(FC) -E -P -x c

FFLAGS = -O2 -ip
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
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
#
#------ Clear this if you do not want netCDF
NETCDF_INCFLAGS= -I$(NETCDF_ROOT)/include
FPPFLAGS_CDF=-DCDF
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdff
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI

LIBS_CPLUS=-lstdc++ 
SCALAPACK_LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

LIBS= $(ELSI_LIB) $(SCALAPACK_LIBS) $(NETCDF_LIBS) $(LIBS_CPLUS)
INCFLAGS += $(ELSI_INCFLAGS)

SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) $(FPPFLAGS_ELSI) -DMPI_TIMING -DF2003 
#
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
