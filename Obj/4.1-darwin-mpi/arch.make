.SUFFIXES:
.SUFFIXES: .f .F .o .c .a .f90 .F90

SIESTA_ARCH = macos-high-sierra

CC = gcc

FPP = $(FC) -E -P -x c

FC = mpif90 -DMPI

FC_SERIAL = gfortran

FFLAGS = -O2 -fPIC -ftree-vectorize -march=native -DFC_HAVE_FLUSH -DFC_HAVE_ABORT

AR = ar
ARFLAGS_EXTRA =

RANLIB = ranlib

SYS = nag

SP_KIND = 4
DP_KIND = 8
KINDS = $(SP_KIND) $(DP_KIND)

DEFS_PREFIX =

LDFLAGS = # -static
FPPFLAGS= -DFC_HAVE_FLUSH -DFC_HAVE_ABORT

FCFLAGS_fixed_f =  $(FFLAGS)
FCFLAGS_free_f90 = $(FFLAGS)
FPPFLAGS_fixed_F = $(FPPFLAGS)
FPPFLAGS_free_F90 = $(FPPFLAGS)

BLAS_LIBS = -lblas
LAPACK_LIBS = -llapack
SCALAPACK_LIBS = -lscalapack -L/usr/local/Cellar/scalapack/2.1.0/lib/

COMP_LIBS =

# By default NetCDF is not used, but it may become mandatory in
# later versions of SIESTA.
#NETCDF_ROOT = /opt/netcdf/4.4.0
#NETCDF_INCFLAGS = -I$(NETCDF_ROOT)/include
#NETCDF_LIBS = -L$(NETCDF_ROOT)/lib -lnetcdff -lnetcdf

MPI_INTERFACE = libmpi_f90.a
MPI_INCLUDE = .

FPPFLAGS = $(DEFS_PREFIX)-DFC_HAVE_ABORT -DHAVE_MPI -DHAVE_SCALAPACK

LIBS = $(NETCDF_LIBS) $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(MPI_LIBS) $(COMP_LIBS)

FFLAGS_DEBUG = -g -O0   # your appropriate flags here...

atom.o: atom.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(CPPFLAGS) $< 
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

