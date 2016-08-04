#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90

PEXSI_VERSION=0.9.0
PEXSI_LIB_DIR=/gpfs/projects/bsc21/bsc21308/SIESTA/build/pexsi_v0.9.0
FPPFLAGS_PEXSI=-DSIESTA__PEXSI
PEXSI_INCFLAGS=-I$(PEXSI_LIB_DIR)/fortran

SIESTA_ARCH=MareNostrum3-intel-openmpi-pexsi$(PEXSI_VERSION)

FPP=
FPP_OUTPUT= 
FC=mpif90
RANLIB=ranlib

SYS=nag

SP_KIND=4
DP_KIND=8
KINDS=$(SP_KIND) $(DP_KIND)

FFLAGS= -O3
FPPFLAGS= -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DMPI_TIMING $(FPPFLAGS_PEXSI)
LDFLAGS=-Vaxlib

ARFLAGS_EXTRA=

FCFLAGS_fixed_f=
FCFLAGS_free_f90=
FPPFLAGS_fixed_F=
FPPFLAGS_free_F90=

MKLROOT=/apps/INTEL/mkl

BLACS_LIBS=-Wl,-rpath,/apps/INTEL/mkl/lib/intel64/ -L/apps/INTEL/mkl/lib/intel64/ -lmkl_blacs_openmpi_lp64
SCALAPACK_LIBS=-Wl,-rpath,/apps/INTEL/mkl/lib/intel64/ -L/apps/INTEL/mkl/lib/intel64/ -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential


#DSUPERLU_DIR  = /gpfs/projects/bsc21/bsc21308/libs/SuperLU_DIST_3.3/
DSUPERLU_DIR  = /apps/SUPERLU/3.3/
SCOTCH_DIR    = /apps/SCOTCH/6.0.0
PARMETIS_DIR  = /apps/PARMETIS/4.0.3
METIS_DIR     = /apps/METIS/5.0.2

#EXTRAE_DIR=/apps/CEPBATOOLS/extrae/latest/default/64

METIS_LIB        = ${METIS_DIR}/lib/libmetis.a
#PARMETISLIB      = ${PARMETIS_DIR}/libparmetis.a
DSUPERLU_LIB     = ${DSUPERLU_DIR}/lib/libsuperlu_dist_3.3.a
PEXSI_LIB        = ${PEXSI_LIB_DIR}/src/libpexsi_release_mn_v$(PEXSI_VERSION).a
SCOTCH_LIB       = -L${SCOTCH_DIR}/lib -lptscotchparmetis -lptscotch -lscotch -lptscotcherr
#EXTRAE_LIB       = -L${EXTRAE_DIR}/lib -lmpitrace

#NETCDF_LIBS=
#NETCDF_INTERFACE=

#LIBS             = $(SCALAPACK_LIBS) $(BLACS_LIBS) ${PEXSI_LIB} ${DSUPERLU_LIB} ${SELINV_LIB} ${PARMETISLIB} ${METIS_LIB} -lstdc++
LIBS             = ${EXTRAE_LIB} $(SCALAPACK_LIBS) $(BLACS_LIBS) ${PEXSI_LIB} ${DSUPERLU_LIB} ${SELINV_LIB} ${SCOTCH_LIB} ${METIS_LIB} -lstdc++

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
