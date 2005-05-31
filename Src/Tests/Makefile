#
# Makefile for SIESTA tests
#
# Location of siesta executable (seen from two levels down)
#
SIESTA=../../../siesta
#SIESTA="/opt/mpich-1.2.6-nag/bin/mpirun -np 2 ../../../siesta"
#
#
###include ../arch.make
#
tests= h2o h2o_dos h2o_orderN batio3 sih si64 fe floating bessel 
#
all:  ###   cmpsig   ## Optional: better not activate as of 1.3
	@for i in ${tests} ; do \
		echo ">>>> Running $$i test..." ; \
		(SIESTA=$(SIESTA) sh do.sh $$i) ; \
	done
#
##cmpsig: cmpsig.f
##	$(FC) $(FFLAGS) -o cmpsig cmpsig.f
#
xmlcheck: all
	@for i in ${tests} ; do \
	  echo "xmllint check $$i.xml ..."; \
	  xmllint $$i/$$i.xml > /dev/null ; \
	done
#
clean:
##	rm -f cmpsig
	@for i in ${tests} ; do \
		echo ">>>> Cleaning $$i test..." ; \
		(cd $$i ; rm -rf work ; rm -f *.out *.xml completed) ; \
	done






