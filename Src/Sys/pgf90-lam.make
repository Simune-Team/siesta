SIESTA_ARCH=hreidar-lam
#
# For hreidar opteron cluster at ETHZ
#
FC=pgf90 -Msecond_underscore
FC_ASIS=$(FC)
#
FFLAGS= -fast
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=dc_lapack.a
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
DEFS_CDF=            #  -DCDF
#
MPI_LIBS= -L/usr/local/apli/lam/lib -llammpio -llamf77mpi -lmpi -llam -lpthread
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/apli/lam/include
DEFS_MPI=-DMPI
#
# There are (were?) some problems with command-line processing compatibility
# that forced the extraction of "pgi.aux" and "pgiarg" as independent 
# libraries (details unfortunately lost)
#
LIBS= -L$(HOME)/lib -lblacs-all-hreidar.lam-0 -lscalapack-hreidar-lam \
      -lblacs-all-hreidar.lam-0 \
      -L/usr/lib64  -llapack -lblas  $(NETCDF_LIBS) $(MPI_LIBS) -lg2c
SYS=cpu_time
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
#
# Important (at least for V5.0-1 of the pgf90 compiler...)
# Compile atom.f without optimization.
#
atom.o:
	$(FC) -c $(FFLAGS_DEBUG) atom.f
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#
