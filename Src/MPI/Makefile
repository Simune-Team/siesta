#
.SUFFIXES: .o .F90 .f90 .F .f .a 
#
include ../arch.make
#
default: module_built
#
TEMPLATES= mpi__type_s.f90 mpi__type_sv.f90 mpi__type_v.f90 mpi__type_vs.f90

OBJS=Interfaces.o
#
SOURCES=$(OBJS:.o=.f90)
#
generated: $(TEMPLATES) 
	@make clean
	@echo
	@echo "PARALLEL VERSION:   **** MPI  WARNING ****"
	@echo "We are assuming that the 'kind' numbers for"
	@echo "single and double precision real are 4 and 8, respectively."
	@echo "This is the case for all compilers we have dealt with so far,"
	@echo "but you should change it if the assumption does not hold for"
	@echo "your compiler!!"
	@echo
	@echo "Hit ^C to abort and change the defaults in MPI/Makefile..."
	@sleep 2
#
	@sh generate.sh "4 8"
	@touch generated
#
$(OBJS): generated  mpi__include.o
	$(FC) -c $(FREE_F90) Interfaces.f90
#
mpi.o:  $(OBJS) mpi__include.o
	$(FC) -c mpi.F
#
module_built: mpi.o $(OBJS) mpi__include.o 
	@ar q libmpi_f90.a mpi.o mpi__include.o $(OBJS)
	-$(RANLIB) libmpi_f90.a
	@mv libmpi_f90.a ..
	@mv *.mod ..
	@touch module_built
#
mpi__include.o: mpi__include.f90
	$(FC) -c $(FREE_F90) -I$(MPI_INCLUDE) mpi__include.f90
#
clean:
	@rm -f $(SOURCES) $(OBJS) mpi.o mpi__include.o
	@rm -f module_built generated *.uses



