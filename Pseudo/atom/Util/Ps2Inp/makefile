.SUFFIXES:
.SUFFIXES: .F90 .f .F .f90 .o .a 
default: ps2inp
#
FC=gfortran
#
OBJS=die.o periodic_table.o ground_state.o flib_spline.o pseudopotential.o make_ps_input.o
ps2inp: $(OBJS)
	gfortran -o ps2inp $(OBJS)
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
