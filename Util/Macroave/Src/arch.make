#
FC=f90 -n32
FFLAGS= -O3
#
.f.o:
	$(FC) -c $(FFLAGS)   $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#
