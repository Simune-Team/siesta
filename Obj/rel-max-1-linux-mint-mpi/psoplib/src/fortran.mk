#
# Fortran macros 
#
FC=gfortran
FFLAGS= -O0 -g -fbacktrace 
LDFLAGS=     
#
AR=ar
RANLIB=echo
#
DEFS=
#
.F.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(INCFLAGS) $(INC_SEARCH) $(FFLAGS)   $<








