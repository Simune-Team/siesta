#
.SUFFIXES: .o .f90 .F .a 
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
	@sh generate.sh "4 8"
	@touch generated
#
$(OBJS): generated  mpi__include.o
	$(FC) -c Interfaces.f90
#
mpi.o:  $(OBJS) mpi__include.o
	$(FC) -c mpi.F
#
module_built: mpi.o $(OBJS) mpi__include.o 
	@ar q libmpi_f90.a mpi.o mpi__include.o $(OBJS)
	@mv libmpi_f90.a ..
	@mv *.mod ..
	@touch module_built
#
mpi__include.o: mpi__include.f90
	$(FC) -c -I$(MPI_INCLUDE) mpi__include.f90
#
clean:
	@rm -f $(SOURCES) $(OBJS) mpi.o mpi__include.o
	@rm -f module_built generated *.uses



