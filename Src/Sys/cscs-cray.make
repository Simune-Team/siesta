SIESTA_ARCH=cscs-cray
#
# For Cray XT-3 at CSCS
#
FC=ftn -target=catamount
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
MPI_LIBS=
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
# There are (were?) some problems with command-line processing compatibility
# that forced the extraction of "pgi.aux" and "pgiarg" as independent 
# libraries (details unfortunately lost)
#
LIBS= 
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
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
