# Edison machine at NERSC
#
PEXSI_VERSION=v0.9.0
#–
SIESTA_ARCH=edison-$(PEXSI_VERSION)
#
FC=ftn
FC_ASIS=$(FC)
#
#FFLAGS= -O3
FFLAGS= -O2
#FFLAGS = -g -O0
FFLAGS_DEBUG= -g -O0
RANLIB=echo

#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
FPPFLAGS_CDF=  -DCDF
#
DUMMY_FOX=--enable-dummy
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=$(MPICH_DIR)/include
KINDS="4 8"
FPPFLAGS_MPI=-DMPI
#
PEXSI_INCFLAGS = -I/project/projectdirs/m1027/PEXSI/libpexsi_interfaces/$(PEXSI_VERSION)
FPPFLAGS_PEXSI=-DPEXSI_$(PEXSI_VERSION)
#
# Extended interface
PEXSI_DIR     = /project/projectdirs/m1027/PEXSI/libpexsi_edison
DSUPERLU_DIR  = /project/projectdirs/m1027/PEXSI/libpexsi_edison
PARMETIS_DIR  = /project/projectdirs/m1027/PEXSI/libpexsi_edison
SCOTCH_DIR    = /project/projectdirs/m1027/PEXSI/libpexsi_edison

METIS_LIB        = ${PARMETIS_DIR}/libmetis.a
SCOTCH_LIB       = -L${SCOTCH_DIR} -lptscotchparmetis -lptscotch -lptscotcherr -lscotch
SELINV_LIB       = ${PEXSI_DIR}/libselinv.a
#DSUPERLU_LIB     = ${DSUPERLU_DIR}/
DSUPERLU_LIB     = ${DSUPERLU_DIR}/libsuperlu_dist_3.3_O2_NOPRNT.a
PEXSI_LIB        = ${PEXSI_DIR}/libpexsi_$(PEXSI_VERSION).a
MKL_LIB          = -mkl=cluster
CPP_LIB          = -lstdc++
LIBS             = ${PEXSI_LIB} ${DSUPERLU_LIB} ${SELINV_LIB} ${SCOTCH_LIB} ${METIS_LIB}\
                   ${MKL_LIB}  ${IPM} ${CPP_LIB}
#
LDFLAGS  = -no-ipo
#
SYS=cpu_time
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) -DMPI_TIMING -DF2003
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
