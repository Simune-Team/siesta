#
# Makefile for SIESTA tests
#
# Location of siesta executable (seen from two levels down)
#
SIESTA=../../../siesta
#
#
tests= h2o h2o_dos h2o_orderN floating bessel si mgco3 \
       var_cell constant_volume batio3 fe fe_broyden sih si64
#
all: 
	@for i in ${tests} ; do \
		(cd $$i ; make SIESTA=$(SIESTA) completed) ; \
	done
#
xmlcheck: all
	@for i in ${tests} ; do \
		(cd $$i ; make xmlcheck) ; \
	done
#
clean:
	@for i in ${tests} ; do \
		(cd $$i ; make clean) ; \
	done






