#
.SUFFIXES: .o .F90 .f90 .F .f .a 
#
default: module_built 
include ../arch.make
#
#
TEMPLATES= mpi__type_s.f90 mpi__type_sv.f90 mpi__type_v.f90 mpi__type_vs.f90

OBJS=Interfaces.o
#
SOURCES=$(OBJS:.o=.f90)
#
generated: $(TEMPLATES)
	@$(MAKE) clean
	@$(MAKE) kind_explorer
	@echo "The kind numbers for single and double precision reals follow"
	@sh generate.sh 
	@sleep 2
	@touch generated
#
kind_explorer: kind_explorer.o
	$(FC) -o kind_explorer kind_explorer.o

$(OBJS): generated  mpi__include.o
	$(FC) $(FFLAGS_MPI)  $(DEFS) -c $(FREE_F90) Interfaces.f90
#
mpi.o:  $(OBJS) mpi__include.o
#
module_built: mpi.o $(OBJS) 
	@ar q libmpi_f90.a mpi.o mpi__include.o $(OBJS)
	-$(RANLIB) libmpi_f90.a
	@mv libmpi_f90.a ..
	@cp *.mod ..
	@touch module_built
#
mpi__include.o: mpi__include.f90
	$(FC) $(DEFS) -c $(FREE_F90) $(FFLAGS_MPI) -I$(MPI_INCLUDE) mpi__include.f90
#
pi3:    pi3.o
	$(FC) -o pi3 $(LDFLAGS) pi3.o ../$(MPI_INTERFACE) $(MPI_LIBS) 
#
blacs_prb:    blacs_prb.o
	$(FC) -o blacs_prb $(LDFLAGS) blacs_prb.o ../$(MPI_INTERFACE) $(LIBS) 
#
pblas_prb:    pblas_prb.o
	$(FC) -o pblas_prb $(LDFLAGS) pblas_prb.o ../$(MPI_INTERFACE) $(LIBS) 
#
clean:
	@rm -f $(SOURCES) $(OBJS) mpi.o mpi__include.o *.mod
	@rm -f module_built generated *.uses kind_explorer
	@rm -f pi3.o pi3




