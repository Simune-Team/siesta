SIESTA_ARCH=nag-darwin-nolibs
#
# Issues addressed: 
#
#      - Case-insensitive filesystem
#      - .F and .F90 suffixes not recognized by f95 on Mac (possibly related
#        to the first point.
#      - CPUtime now called in the f95 fashion (file nag.f)
#      - -dusty option to allow EQUIVALENCE statement in blas.f
#      - -dcfuns option to allow non-standard intrinsics such as dconjg
# 
# Other issues:
#  
#      - Ranlib is necessary for the .a files. As of now it has to be
#        done manually. Later on, one could include a statement such as
#
#                        -$(RANLIB) XXXX.a
#
#        in the makefiles for fdf and Libs. But sometimes this does not
#        work if the file is later moved...
#
#      ONE HAS TO BE EXTRA CAREFUL WHEN DELETING FILES, DUE TO THE
#      CASE-INSENSITIVENESS OF THE FILE SYSTEM ON THE MAC.
#
FC=f95
FC_ASIS=$(FC)
RANLIB=ranlib
#
FFLAGS= -g -dcfuns -dusty
FFLAGS_DEBUG= 
LDFLAGS=
COMP_LIBS=
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS= 
COMP_LIBS=linalg.a
SYS=nag
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
CPP=/usr/local/lib/NAGWare/fpp -P
#
.F.o:
	mv $<  $*.fpp
	$(CPP) -fixed $(DEFS) $*.fpp > $*.f
	$(FC) -c $(FFLAGS)  $(DEFS) $*.f
	@rm -f $*.f
	@mv $*.fpp $*.F
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	mv $<  $*.fpp90
	$(CPP) -free $(DEFS) $*.fpp90 > $*.f90
	$(FC) -c $(FFLAGS) $*.f90
	@rm -f $*.f90
	@mv $*.fpp90 $*.F90
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#









