#
.SUFFIXES: .o .f90 .a 
#
include ../Sys/${SIESTA_SYS}.make
#
default: module_built
#
TEMPLATES= mpi__type_s.f90 mpi__type_sv.f90 mpi__type_v.f90 mpi__type_vs.f90

OBJS=mpi__character_s.o mpi__character_sv.o mpi__character_v.o \
     mpi__character_vs.o mpi__complex_s.o mpi__complex_sv.o \
     mpi__complex_v.o mpi__complex_vs.o mpi__dcomplex_s.o \
     mpi__dcomplex_sv.o mpi__dcomplex_v.o mpi__dcomplex_vs.o \
     mpi__double_s.o mpi__double_sv.o mpi__double_v.o \
     mpi__double_vs.o mpi__integer_s.o \
     mpi__integer_sv.o mpi__integer_v.o mpi__integer_vs.o \
     mpi__logical_s.o mpi__logical_sv.o mpi__logical_v.o \
     mpi__logical_vs.o mpi__single_s.o mpi__single_sv.o \
     mpi__single_v.o mpi__single_vs.o 
#
SOURCES=$(OBJS:.o=.f90)
#
.generated: $(TEMPLATES)
	@make clean
	@sh generate.sh
	@touch .generated
#
$(OBJS): .generated  mpi__include.o
	$(FC) -c $*.f90
#
mpi.o:  $(OBJS) mpi__include.o
	$(FC) -c mpi.f90
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
	@rm -f $(SOURCES) *.o ../libmpi_f90.a 
	@rm -f ../mpi*.mod 
	@rm -f module_built .generated



