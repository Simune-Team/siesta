#
# Makefile for SIESTA tests
#
# Location of siesta executable (seen from two levels down)
#
SIESTA=../../../siesta
#SIESTA="mpirun -np 4 ../../../siesta"
#
#
include ../arch.make
#
tests= h2o  batio3 sih si64 fe # h2o_orderN 
#
all: cmpsig
	@for i in ${tests} ; do \
		echo ">>>> Running $$i test..." ; \
		(SIESTA=$(SIESTA) sh do.sh $$i) ; \
	done
#
cmpsig: cmpsig.f
	$(FC) $(FFLAGS) -o cmpsig cmpsig.f
#
clean:
	rm -f cmpsig
	@for i in ${tests} ; do \
		echo ">>>> Cleaning $$i test..." ; \
		(cd $$i ; rm -rf work ; rm -f *.out completed) ; \
	done
