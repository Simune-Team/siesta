SIESTA_ARCH=cscs-cray-mpi
#
# For Cray XT-3 at CSCS
#
FC=ftn -target=catamount
FPP=linux-pgf90 -F
FC_ASIS=$(FC)
#
FFLAGS= -g -fastsse
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=dc_lapack.a
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
DEFS_CDF=            #  -DCDF
#
KINDS=4 8
MPI_LIBS=
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.
DEFS_MPI=-DMPI
#
LIBS= 
SYS=nag
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
#
# Important 
# Compile atom.f and electrostatic.f without optimization.
#
atom.o:
	$(FC) -c $(FFLAGS_DEBUG) atom.f
electrostatic.o:
	$(FC) -c $(FFLAGS_DEBUG) electrostatic.f
#
.F.o:
	$(FPP) $(DEFS) $<  ; mv $*.f aux_$*.f
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) -o $*.o aux_$*.f
	rm -f aux_$*.f
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FPP) $(DEFS) $<  ; mv $*.f aux_$*.f90
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) -o $*.o aux_$*.f90
	rm -f aux_$*.f90
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
