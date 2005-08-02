#
# Makefile for SIESTA tests
#
# Location of siesta executable (seen from two levels down)
#
SIESTA=../../../siesta
#
#
tests= h2o h2o_dos h2o_orderN batio3 fe fe_broyden sih si64  floating bessel 
#
all: 
	@for i in ${tests} ; do \
		echo ">>>> Running $$i test..." ; \
		(SIESTA=$(SIESTA) sh do.sh $$i) ; \
	done
#
xmlcheck: all
	@for i in ${tests} ; do \
	  echo "xmllint check $$i.xml ..."; \
	  xmllint $$i/$$i.xml > /dev/null ; \
	done
#
clean:
	@for i in ${tests} ; do \
		echo ">>>> Cleaning $$i test..." ; \
		(cd $$i ; rm -rf work ; rm -f *.out *.xml completed) ; \
	done






