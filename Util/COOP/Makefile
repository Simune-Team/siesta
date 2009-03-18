.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o
#
# Choose compiler here.
#
#FC=g95 -g -O0 -fbounds-check
#FC=ifort -g -O1 -CB -ftrapuv
#FC=ifort -openmp -O3
#FC=ifort  -O3
FC=gfortran -O3 -DGFORTRAN
#FC=gfortran  -openmp -O3 -DGFORTRAN
#
MPROP_OBJS=m_getopts.o f2kcli.o subs.o precision.o units.o \
           io.o sys.o read_curves.o io_hs.o orbital_set.o main_vars.o
#
default: mprop
#
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
#
mprop: $(MPROP_OBJS) mprop.o
	$(FC) -o $@ $(MPROP_OBJS) mprop.o
#
clean:
	rm -f *.o *.mod mprop 

.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#


# DO NOT DELETE THIS LINE - used by make depend
io.o: sys.o
io_hs.o: main_vars.o precision.o
m_getopts.o: f2kcli.o
main_vars.o: m_getopts.o precision.o subs.o units.o
mprop.o: io_hs.o main_vars.o orbital_set.o read_curves.o
orbital_set.o: main_vars.o subs.o
read_curves.o: orbital_set.o precision.o
subs.o: precision.o
units.o: precision.o
