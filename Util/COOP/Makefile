.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o
#
# Choose compiler here.
# If you have netCDF installed, do something similar to this:
#NETCDF_ROOT=$(HOME)/lib/netcdf-3.6.2-gfortran
#INCFLAGS=-I$(NETCDF_ROOT)/include
#
#DEFS=-DGFORTRAN -DCDF 
#
#NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
#
#FC=gfortran -O3 $(DEFS) $(INCFLAGS)
#
# ---------------------
FC=gfortran -O3 -DGFORTRAN
# ---------------------
#
MPROP_OBJS=m_getopts.o f2kcli.o subs.o precision.o units.o \
           io.o sys.o read_curves.o io_hs.o orbital_set.o main_vars.o
DM_OBJS=m_getopts.o f2kcli.o subs.o precision.o units.o \
           io.o sys.o iodm_netcdf.o io_hs.o  main_vars.o write_dm.o
#
default: mprop
#
dep: 
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90
#
mprop: $(MPROP_OBJS) mprop.o
	$(FC) -o $@ $(MPROP_OBJS) mprop.o
dm_creator: $(DM_OBJS) dm_creator.o
	$(FC) -o $@ $(DM_OBJS) dm_creator.o $(NETCDF_LIBS)
#
clean:
	rm -f *.o *.mod mprop  dm_creator

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
dm_creator.o: io_hs.o iodm_netcdf.o main_vars.o
io.o: sys.o
io_hs.o: main_vars.o precision.o
iodm.o: precision.o sys.o
iodm_netcdf.o: precision.o sys.o
m_getopts.o: f2kcli.o
main_vars.o: m_getopts.o precision.o subs.o units.o
mprop.o: io_hs.o main_vars.o orbital_set.o read_curves.o
old_dm.o: io_hs.o iodm_netcdf.o main_vars.o orbital_set.o
orbital_set.o: main_vars.o subs.o
read_curves.o: orbital_set.o precision.o
subs.o: precision.o
units.o: precision.o
write_dm.o: precision.o sys.o
