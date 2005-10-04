.SUFFIXES: .o .F90 .f90 .F .f .a 
#
default: module_built 
include ../arch.make
#
#
TEMPLATES= mpi__type_s.f90 mpi__type_sv.f90 mpi__type_v.f90 mpi__type_vs.f90
#
INCFLAGS=-I$(MPI_INCLUDE)
#
Interfaces.f90 V_S.uses VS.uses: $(TEMPLATES)
	if [ -z "$(KINDS)" ] ; then  $(MAKE) kind_explorer ; fi
	@echo "The kind numbers for single and double precision reals follow"
	@sh generate.sh "$(KINDS)"
#
kind_explorer: kind_explorer.o
	$(FC) -o kind_explorer $(LDFLAGS) kind_explorer.o
mpi.o: mpi__include.o Interfaces.o V_S.uses VS.uses
#
module_built: mpi.o mpi__include.o Interfaces.o
	@ar $(ARFLAGS_EXTRA) cru libmpi_f90.a mpi.o mpi__include.o Interfaces.o
	-$(RANLIB) libmpi_f90.a
	@mv libmpi_f90.a ..
	@cp *.mod ..
	@touch module_built
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
	@rm -f Interfaces.f90 *.o *.mod
	@rm -f module_built *.uses kind_explorer
	@rm -f pi3 blacs_prb pblas_prb




