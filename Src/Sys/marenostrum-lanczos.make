BIT=64
SIESTA_ARCH=XLF $(BIT)bits PARALLEL MPI
ROGELI_HOME=/home/bsc21/bsc21017
#
FC= mpif90   #xlf90_r
CC= mpicc
#
FFLAGS_DEBUG= -g
#FFLAGS=-O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q$(BIT)
#FFLAGS=-O3 -qstrict -qarch=com -q$(BIT)
#FFLAGS=-g -O0 -qfullpath -q$(BIT) -qarch=com

FFLAGS=-O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q$(BIT)
CFLAGS=-O3 -qtune=ppc970 -qarch=ppc970 -q$(BIT)

#FFLAGS_parse=-qsuffix=f=f -qfree #-qfixed
LDFLAGS= -q$(BIT)
#LDFLAGS= -q$(BIT) -g
COMP_LIBS=
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=    libmpi_f90.a
MPI_INCLUDE=.
MPI_LIBS=         #-lblacsgm
DEFS_MPI=         -WF,-DMPI



#BLAS = -L/gpfs/apps/LAPACK/lib$(BIT) -llapack \
#	     -L/gpfs/apps/SCALAPACK/lib$(BIT) -lblas


BLAS= -L/gpfs/apps/SCALAPACK/lib$(BIT) -lscalapack \
        /gpfs/apps/SCALAPACK/lib$(BIT)/blacsF77init_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib$(BIT)/blacsCinit_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib$(BIT)/blacs_MPI-PPC-0.a \
      -L/gpfs/apps/LAPACK/lib$(BIT) -llapack \
      -lessl
#      -L/gpfs/apps/SCALAPACK/lib$(BIT) -lblas

METIS= -L/gpfs/apps/METIS/4.0/$(BIT)/ -lmetis

#TRACE_DIR=/gpfs/apps/CEPBATOOLS/mpitrace-mx/$(BIT)/lib/
#MPITRACER=-L$(TRACE_DIR) -lmpitracef -lxml2
#MPITRACER=$(TRACE_DIR)/libmpitracef.a -lxml2 /opt/osshpc/mpich-mx/64/lib/libmpich.a
#DEFS_TRACE= -WF,-DMPI_TRACE
#PAPI=-L/gpfs/apps/PAPI/papi-3.2.1/$(BIT)/lib/ -lpapi

#MPITRACE_HOME = /home/bsc41/bsc41127/Tools/mpitrace_mrnet_18_04_2008/Package
#PAPI_HOME = /gpfs/apps/PAPI/3.2.1-970mp/64
#MRNET_HOME = /home/bsc41/bsc41127/apps/mrnet_last
#MPITRACER = $(MPITRACE_HOME)/lib/libmpitracef.a -L$(MRNET_HOME)/lib -lmrnet -lxplat -lxml2 -ldl -lpthread -lstdc++
#PAPI=$(PAPI_HOME)/lib/libpapi.a $(PAPI_HOME)/lib/libperfctr.a

PARMETIS=$(ROGELI_HOME)/ParMetis-3.1.1/libparmetis.a \
         $(ROGELI_HOME)/ParMetis-3.1.1/libmetis.a

SLU_HOME=$(ROGELI_HOME)/SuperLU_DIST_2.3
SLU=$(SLU_HOME)/lib/libsuperlu_dist_2.3.a
INCLUDE=-I$(SLU_HOME)/SRC

ARPACK_LIB=$(ROGELI_HOME)/ARPACK/parpack_MPI-PPC.a $(ROGELI_HOME)/ARPACK/libarpack_PPC.a

OMP_LIB=
OMP_DEF=

SYS=xlf

LIBS= $(SLU) $(PARMETIS) $(METIS) $(ARPACK_LIB) $(BLAS) $(PAPI) \
      $(OMP_LIB) $(MPITRACER) $(-L/opt/osshpc/mx/lib$(BIT)/

DEFS= $(DEFS_MPI) $(DEFS_CDF) $(DEFS_TRACE) $(OMP_LIB) -WF,-DDEBUG -WF,-DXLF

FREE_F90=-qsuffix=f=f90
#
.c.o:
	$(CC) -c -DNoChange $(CFLAGS) $(INCLUDE) $(CDEFS) $<

.F90.o:
	$(FC) -qsuffix=cpp=F90 -c $(INCFLAGS) $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -qsuffix=f=f90 -c $(INCFLAGS) $(FFLAGS)   $<
.F.o:
	$(FC) -qsuffix=cpp=F -c $(INCFLAGS) -qfixed $(FFLAGS) $(DEFS) $<
#.f.o:
#	$(FC) -qsuffix=f=f -qfixed -c $(INCFLAGS) $(FFLAGS)   $<
#
.f.o:
	$(FC) -qsuffix=cpp=f -c $(INCFLAGS) -qfixed $(FFLAGS) $(DEFS) $<

