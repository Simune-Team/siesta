SIESTA_ARCH=osfdxml
#
FC=f90
FFLAGS= -fast -tune host
FFLAGS_DEBUG= -g3 -fast -tune host
LIBS= -ldxml
SYS=bsd
DEFS=
MPILIB=
CPP=/bin/cpp -P
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(CPP) $(DEFS) $< > $*.f90
	$(FC) -c $(FFLAGS) $*.f90
	@rm -f $*.f90
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#
