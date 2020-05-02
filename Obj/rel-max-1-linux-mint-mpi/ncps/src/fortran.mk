#
# Fortran macros 
#
FC=gfortran
FFLAGS= -O0 -g -fbacktrace 
LDFLAGS=     
#
PSML_ROOT=$(HOME)/code/SIESTA/S2.6/GIT/psml/src
PSML_INCFLAGS=-I$(PSML_ROOT)
PSML_LIBS=$(PSML_ROOT)/libpsml.a
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








