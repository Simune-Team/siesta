#
#  Makefile for FDF library and FDF example.
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
LOCAL_INCFLAGS=-I../
include ../arch.make
INCFLAGS:=$(INCFLAGS) $(LOCAL_INCFLAGS)
#
# Library module libfdf.a
#
module: libfdf.a
libfdf.a: iso_fortran_env.o fdf.o io_fdf.o parse.o utils.o prec.o 
	$(AR) $(ARFLAGS_EXTRA) cru $@ $^
	-$(RANLIB) $@
	cp libfdf.a ..
	@cp -p *.*d ..     

fdf.o: io_fdf.o parse.o prec.o utils.o
io_fdf.o: utils.o prec.o
parse.o: utils.o prec.o
utils.o: prec.o

##%.o: %.F90
##	$(FC) $(FFLAGS) -c $< $(DEFS)

#
# Tests for FDF library (Sequential, MPI, IO)
#
sample: sample.F90 libfdf.a
	$(FC) $(LDFLAGS) -o $@ $^ $(MPIFLAGS)

#
# Cleaning
#
clean:
	rm -f *.o *.mod libfdf.a
	rm -f sample sample_mpi sample_cluster sample_blocking
	rm -f fdf.debug* sample.out* sample.fdf.*
